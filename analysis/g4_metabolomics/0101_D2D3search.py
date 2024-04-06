## Holder for D2D3search.py
import sys
sys.path.append("./")
import numpy as np
import pandas as pd
from corems.mass_spectra.input.mzml import MZMLSpectraParser
from corems.mass_spectra.output.export import Lipidomics_Export
from pathlib import Path
from scipy.spatial import KDTree
from scipy import sparse


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
    myLCMSobj.parameters.lc_ms.ph_inten_min = myLCMSobj._ms_unprocessed[1].intensity.max()*thresh
    myLCMSobj.parameters.lc_ms.ph_persis_min = myLCMSobj._ms_unprocessed[1].intensity.max()*thresh
    myLCMSobj.parameters.lc_ms.ph_smooth_it = 0
    myLCMSobj.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel = 5*10**-6
    myLCMSobj.parameters.lc_ms.mass_feature_cluster_rt_tolerance = 1

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
    myLCMSobj.find_mass_features(verbose = verbose) 
    myLCMSobj.integrate_mass_features(drop_if_fail = True)

def find_H_isos(myLCMSobj, D):
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

    mz_diff = 1.006276745946*D # D2 mass difference
    tol = [mf_df['mz'].median()*myLCMSobj.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel , myLCMSobj.parameters.lc_ms.mass_feature_cluster_rt_tolerance*0.5]  # mz, in relative; scan_time in minutes

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
    pairs = np.stack((distances.row, distances.col), axis=1)  # 1H to 2H pairs

    # Turn pairs (which are index of mf_df) into mf_id and then into two dataframes to join to mf_df
    pairs_mf = pairs.copy()
    pairs_mf[:,0] = mf_df.iloc[pairs[:,0]].mf_id.values
    pairs_mf[:,1] = mf_df.iloc[pairs[:,1]].mf_id.values

    monos = np.setdiff1d(np.unique(pairs_mf[:, 0]), np.unique(pairs_mf[:, 1])) 
    for mono in monos:
        myLCMSobj.mass_features[mono].monoisotopic_mf_id = mono
    pairs_iso_df = pd.DataFrame(pairs_mf, columns=['parent', 'child'])
    while not pairs_iso_df.empty:
        pairs_iso_df = pairs_iso_df.set_index('parent', drop=False)
        m1_isos = pairs_iso_df.loc[monos, 'child'].unique()
        for iso in m1_isos:
            # Set monoisotopic_mf_id and isotopologue_type for isotopologues
            parent = pairs_mf[pairs_mf[:, 1] == iso, 0]
            if len(parent) > 1:
                # Choose the parent that is closest in time to the isotopologue
                parent_time = [myLCMSobj.mass_features[p].scan_time for p in parent]
                time_diff = [np.abs(myLCMSobj.mass_features[iso].scan_time - x) for x in parent_time]
                parent = parent[np.argmin(time_diff)]
            else:
                parent = parent[0]
            myLCMSobj.mass_features[iso].monoisotopic_mf_id = myLCMSobj.mass_features[parent].monoisotopic_mf_id
            if myLCMSobj.mass_features[iso].monoisotopic_mf_id is not None:
                mass_diff = myLCMSobj.mass_features[iso].mz - myLCMSobj.mass_features[myLCMSobj.mass_features[iso].monoisotopic_mf_id].mz
                myLCMSobj.mass_features[iso].isotopologue_type = "2H"+ str(int(round(mass_diff, 0)))

        # Drop the mono and iso from the pairs_iso_df
        pairs_iso_df = pairs_iso_df.drop(index = monos, errors = 'ignore') #Drop pairs where the parent is a child that is a child of a root
        pairs_iso_df = pairs_iso_df.set_index('child', drop=False)
        pairs_iso_df = pairs_iso_df.drop(index = m1_isos, errors = 'ignore')

        if not pairs_iso_df.empty:
            # Get new monos, recognizing that these are just 13C isotopologues that are connected to other 13C isotopologues to repeat the process
            monos = np.setdiff1d(np.unique(pairs_iso_df.parent), np.unique(pairs_iso_df.child))

if __name__ == '__main__':

    # Run positive mode 
    file_dir = Path("data/raw/metabolomics/G4/D3_Homarine_Fate_Inc/mzML/positive")
    out_dir = Path("data/raw/metabolomics/G4/D3_Homarine_Fate_Inc/d3_search_results/positive")
    out_dir.mkdir(parents=True, exist_ok=True)
    files_list = list(file_dir.glob("*.mzML"))
    out_paths_list = [out_dir / f.stem for f in files_list]
    for file_in, file_out in list(zip(files_list, out_paths_list)):
        print(f"Processing {file_in}")
        myLCMSobj = instantiate_lcms_obj(file_in, verbose = True)
        set_params_on_lcms_obj(myLCMSobj, thresh = 0.0001)
        signal_processing_lcms(myLCMSobj, verbose = True)
        find_H_isos(myLCMSobj, D = 2)
        find_H_isos(myLCMSobj, D = 3)
        mf_df = myLCMSobj.mass_features_to_df()
        # drop all mass features that have no monoisotopic_mf_id 
        mf_df = mf_df[mf_df['monoisotopic_mf_id'].notnull()]
        mf_df.to_csv(file_out.with_suffix('.csv'))

    # Run negative mode 
    file_dir = Path("data/raw/metabolomics/G4/D3_Homarine_Fate_Inc/mzML/negative")
    out_dir = Path("data/raw/metabolomics/G4/D3_Homarine_Fate_Inc/d3_search_results/negative")
    out_dir.mkdir(parents=True, exist_ok=True)
    files_list = list(file_dir.glob("*.mzML"))
    out_paths_list = [out_dir / f.stem for f in files_list]
    for file_in, file_out in list(zip(files_list, out_paths_list)):
        print(f"Processing {file_in}")
        myLCMSobj = instantiate_lcms_obj(file_in, verbose = True)
        set_params_on_lcms_obj(myLCMSobj, thresh = 0.0001)
        signal_processing_lcms(myLCMSobj, verbose = True)
        find_H_isos(myLCMSobj, D = 2)
        find_H_isos(myLCMSobj, D = 3)
        mf_df = myLCMSobj.mass_features_to_df()
        # drop all mass features that have no monoisotopic_mf_id 
        mf_df = mf_df[mf_df['monoisotopic_mf_id'].notnull()]
        mf_df.to_csv(file_out.with_suffix('.csv'))

    print("Finished processing files")


