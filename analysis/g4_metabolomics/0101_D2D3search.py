## Holder for D2D3search.py
import sys
sys.path.append("./")
import numpy as np
import pandas as pd
from corems.mass_spectra.input.mzml import MZMLSpectraParser
from pathlib import Path
from scipy.spatial import KDTree
from scipy import sparse
from exporter import LCMSExport


def instantiate_lcms_obj(file_in, verbose):
    """Instantiate a corems LCMS object from a binary file.  Pull in ms1 spectra into dataframe (without storing as MassSpectrum objects to save memory)
    
    Parameters
    ----------
    file_in : str or Path
        Path to binary file
    verbose : bool
        Whether to print verbose output
        
    Returns
    -------
    myLCMSobj : corems LCMS object
        LCMS object with ms1 spectra in dataframe
    """
    # Instantiate parser based on binary file type
    if ".mzML" in str(file_in):
        parser = MZMLSpectraParser(file_in)

    # Instantiate lc-ms data object using parser and pull in ms1 spectra into dataframe (without storing as MassSpectrum objects to save memory)
    myLCMSobj = parser.get_lcms_obj(spectra = "ms1", verbose=verbose)

    return myLCMSobj

def set_params_on_lcms_obj(myLCMSobj, thresh):
    """Set parameters on the LCMS object
    
    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to set parameters on
    thresh : float
        Relative intensity and persistence thresholds (as a fraction of max intensity of a single profile mass)

    Returns
    -------
    None, sets parameters on the LCMS object
    """
    ## persistent homology parameters
    myLCMSobj.parameters.lc_ms.peak_picking_method = "persistent homology"
    myLCMSobj.parameters.lc_ms.ph_inten_min_rel = thresh
    myLCMSobj.parameters.lc_ms.ph_persis_min_rel = thresh
    myLCMSobj.parameters.lc_ms.ph_smooth_it = 0
    myLCMSobj.parameters.lc_ms.ph_smooth_radius_scan = 5
    myLCMSobj.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel = 8*10**-6
    myLCMSobj.parameters.lc_ms.mass_feature_cluster_rt_tolerance = 1.5
    myLCMSobj.parameters.mass_spectrum.noise_threshold_method = "relative_abundance"
    myLCMSobj.parameters.mass_spectrum.noise_threshold_min_relative_abundance = 1
    myLCMSobj.parameters.mass_spectrum.noise_min_mz = 0
    myLCMSobj.parameters.mass_spectrum.noise_max_mz = 2500
    myLCMSobj.parameters.mass_spectrum.min_picking_mz = 0
    myLCMSobj.parameters.mass_spectrum.max_picking_mz = 2500

def signal_processing_lcms(myLCMSobj, verbose):
    """Signal processing for LCMS object.  
    
    This includes peak picking, peak grouping, peak integration, annotation of c13 mass features.

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to process
    verbose : bool
        Whether to print verbose output

    Returns
    -------
    None, processes the LCMS object    
    """
    # Find mass features, cluster, and integrate them.  Then annotate pairs of mass features that are c13 iso pairs.
    myLCMSobj.find_mass_features(verbose=False)
    myLCMSobj.cluster_mass_features(verbose=False)
    myLCMSobj.integrate_mass_features(drop_if_fail=True)

def find_H_isos(myLCMSobj, D = 2):
    """Find D2 pairs in the LCMS object.  
    
    This includes finding pairs of mass features that are D2 isotopes.

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to process
    verbose : bool
        Whether to print verbose output

    Returns
    -------
    None, processes the LCMS object    
    """
    # Data prep fo sparse distance matrix
    dims = ['mz', 'scan_time'] 
    mf_df = myLCMSobj.mass_features_to_df().copy()
    # drop all mass features that have no area (these are likely to be noise)
    mf_df = mf_df[mf_df['area'].notnull()]
    mf_df['mf_id'] = mf_df.index.values
    dims = ['mz', 'scan_time']

    # Sort my ascending mz so we always get the monoisotopic mass first, regardless of the order/intensity of the mass features
    mf_df = mf_df.sort_values(by=['mz']).reset_index(drop=True).copy()
    # Remove monoisotopic_mf_id and isotopolgue_type columns
    mf_df.drop(columns=['monoisotopic_mf_id', 'isotopologue_type'], inplace=True)

    mz_diff = 1.006276745946*D# mass difference if you add two deuteriums
    tol = [mf_df['mz'].median()*myLCMSobj.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel*0.5, myLCMSobj.parameters.lc_ms.mass_feature_cluster_rt_tolerance*0.1]  # mz, in relative; scan_time in minutes

    # Compute inter-feature distances
    distances = None
    for i in range(len(dims)):
        # Construct k-d tree
        values = mf_df[dims[i]].values
        tree = KDTree(values.reshape(-1, 1))

        max_tol = tol[i]
        if dims[i] == 'mz':
            # Maximum absolute tolerance
            max_tol = mz_diff + tol[i]

        # Compute sparse distance matrix
        # the larger the max_tol, the slower this operation is
        sdm = tree.sparse_distance_matrix(tree, max_tol, output_type='coo_matrix')

        # Only consider forward case, exclude diagonal
        sdm = sparse.triu(sdm, k=1)

        if dims[i] == 'mz':
            min_tol = mz_diff - tol[i]
            # Get only the ones that are above the min tol
            idx = sdm.data > min_tol

            # Reconstruct sparse distance matrix
            sdm = sparse.coo_matrix((sdm.data[idx], (sdm.row[idx], sdm.col[idx])),
                                    shape=(len(values), len(values)))

        # Cast as binary matrix
        sdm.data = np.ones_like(sdm.data)

        # Stack distances
        if distances is None:
            distances = sdm
        else:
            distances = distances.multiply(sdm)

    # Extract indices of within-tolerance points
    distances = distances.tocoo()
    pairs = np.stack((distances.row, distances.col), axis=1)  # C12 to C13 pairs

    # Turn pairs (which are index of mf_df) into mf_id and then into two dataframes to join to mf_df
    pairs_mf = pairs.copy()
    pairs_mf[:, 0] = mf_df.iloc[pairs[:, 0]].mf_id.values
    pairs_mf[:, 1] = mf_df.iloc[pairs[:, 1]].mf_id.values

    # Connect monoisotopic masses with isotopologes within mass_features
    pairs_iso_df = pd.DataFrame(pairs_mf, columns=["parent", "child"])
    # drop rows where child number is greater than parent number (indicates more D2 tha original), only if D = 2
    if D == 2:
        pairs_iso_df = pairs_iso_df[pairs_iso_df.child > pairs_iso_df.parent].copy()
    # order by largest parent, then largest daughter
    pairs_iso_df.sort_values(["child", "parent"], inplace=True, ascending=False)
    pairs_iso_df = pairs_iso_df.set_index("parent", drop=False)
    monos = np.setdiff1d(np.unique(pairs_iso_df.parent), np.unique(pairs_iso_df.child))

    m1_isos = pairs_iso_df.loc[monos, "child"].unique()
    for iso in m1_isos:
        # Set monoisotopic_mf_id and isotopologue_type for isotopologues
        parent = pairs_mf[pairs_mf[:, 1] == iso, 0].min()
        # If the feature has no monoisotopic_mf_id, set it and the monoisotopic_mf_id of the isotopologue to the parent
        if myLCMSobj.mass_features[iso].monoisotopic_mf_id is None:
            if myLCMSobj.mass_features[parent].monoisotopic_mf_id is None:
                myLCMSobj.mass_features[parent].monoisotopic_mf_id = parent
            myLCMSobj.mass_features[iso].monoisotopic_mf_id = myLCMSobj.mass_features[parent].monoisotopic_mf_id
            if myLCMSobj.mass_features[iso].monoisotopic_mf_id is not None:
                myLCMSobj.mass_features[iso].isotopologue_type = "2H" + str(D)

def run(files_list, out_files_list):
    for file_in, file_out in list(zip(files_list, out_paths_list)):
        csv_out = out_dir / (file_out.stem + "_isos.csv")
        if csv_out.exists():
            print(f"Skipping {file_out} because it already exists")
            continue
        print(f"Processing {file_in}")
        myLCMSobj = instantiate_lcms_obj(file_in, verbose = True)
        print(f"Instantiated LCMS object from {file_in}")
        set_params_on_lcms_obj(myLCMSobj, thresh = 0.00001)
        try:
            signal_processing_lcms(myLCMSobj, verbose = True)
        except:
            print(f"Failed to process {file_in}")
            continue
        myLCMSobj.find_c13_mass_features()
        
        # Now find D2 isotopologues
        find_H_isos(myLCMSobj, D = 2)
       
        # Now find D3 isotopologues
        find_H_isos(myLCMSobj, D = 3)
        # Write out the mass features to a csv with the D2 and D3 isotopologues
        mf_df = myLCMSobj.mass_features_to_df()
        mf_df.to_csv(out_dir / (file_out.stem + "_isos.csv"))

        # Add the ms1 spectrum to the LCMS object for the mono isotopic mass feature (associated with D2 or D3 mass feature) so we can plot them later
        mono_mf_ids = list(mf_df['monoisotopic_mf_id'].unique())
        # Drop any NAs or None from the list
        mono_mf_ids = [x for x in mono_mf_ids if x is not None]
        for mono_mf_id in mono_mf_ids:
            mono_mf_scan = int(myLCMSobj.mass_features[mono_mf_id].apex_scan)
            myLCMSobj.add_mass_spectra(
                [mono_mf_scan], 
                spectrum_mode="profile",
                use_parser=False)
            myLCMSobj.mass_features[mono_mf_id].mass_spectrum = myLCMSobj._ms[mono_mf_scan]
        
        # Export the lcms object to an hdf5 file
        exporter = LCMSExport(str(file_out), myLCMSobj)
        exporter.to_hdf(overwrite=True)
        print("Exported to hdf5")


if __name__ == '__main__':

    # Run positive mode 
    file_dir = Path("data/raw/metabolomics/G4/D3_Homarine_Fate_Inc/mzML/positive")
    out_dir = Path("data/intermediate/metabolomics/g4/d3d2_search_results/positive")
    out_dir.mkdir(parents=True, exist_ok=True)
    files_list = list(file_dir.glob("*.mzML"))
    files_list = [x for x in files_list if "_Smp_" in x.stem]
    out_paths_list = [out_dir / f.stem for f in files_list]
    run(files_list=files_list, out_files_list=out_paths_list)
    

    # Run negative mode 
    file_dir = Path("data/raw/metabolomics/G4/D3_Homarine_Fate_Inc/mzML/negative")
    out_dir = Path("data/intermediate/metabolomics/g4/d3d2_search_results/negative")
    out_dir.mkdir(parents=True, exist_ok=True)
    files_list = list(file_dir.glob("*.mzML"))
    files_list = [x for x in files_list if "_Smp_" in x.stem]
    out_paths_list = [out_dir / f.stem for f in files_list]
    run(files_list=files_list, out_files_list=out_paths_list)

    print("Finished processing files")


