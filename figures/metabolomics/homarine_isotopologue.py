import math
from rdkit.Chem import rdMolTransforms

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

# Create molecule from SMILES
smiles = "C[N+]1=CC=CC=C1C(=O)[O-]"
mol = Chem.MolFromSmiles(smiles)

# Add explicit hydrogens
mol = Chem.AddHs(mol)

# Generate 2D coordinates
AllChem.Compute2DCoords(mol)


# Rotate molecule so methyl carbon is at top
conf = mol.GetConformer()
# Identify the methyl carbon (C attached to N that only has H neighbors besides N)
methyl_carbon_idx = None
for atom in mol.GetAtoms():
    if atom.GetSymbol() == "N":
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == "C":
                # Check if this carbon only has N and H as neighbors (methyl group)
                neighbor_symbols = [n.GetSymbol() for n in neighbor.GetNeighbors()]
                if neighbor_symbols.count("H") >= 3:  # Methyl group has 3 H's
                    methyl_carbon_idx = neighbor.GetIdx()
                    break
        if methyl_carbon_idx is not None:
            break

if methyl_carbon_idx is None:
    raise ValueError("Could not find methyl carbon atom.")

# Get the center of the molecule
center_x = (
    sum(conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())) / mol.GetNumAtoms()
)
center_y = (
    sum(conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())) / mol.GetNumAtoms()
)

# Calculate angle to rotate methyl carbon to top
methyl_pos = conf.GetAtomPosition(methyl_carbon_idx)
angle = math.atan2(methyl_pos.y - center_y, methyl_pos.x - center_x)
rotation = -angle + math.pi / 2  # Rotate to put methyl carbon at top

# Add 10 degrees (in radians) counterclockwise
rotation += math.radians(-10)

# Apply rotation
for i in range(mol.GetNumAtoms()):
    pos = conf.GetAtomPosition(i)
    x = pos.x - center_x
    y = pos.y - center_y
    new_x = x * math.cos(rotation) - y * math.sin(rotation) + center_x
    new_y = x * math.sin(rotation) + y * math.cos(rotation) + center_y
    conf.SetAtomPosition(i, (new_x, new_y, 0))

# Flip molecule horizontally (mirror x-coordinates around center_x)
for i in range(mol.GetNumAtoms()):
    pos = conf.GetAtomPosition(i)
    flipped_x = 2 * center_x - pos.x
    conf.SetAtomPosition(i, (flipped_x, pos.y, 0))

# Create drawer
drawer = rdMolDraw2D.MolDraw2DCairo(800, 600)

# Identify the methyl carbon (C attached to N that only has H neighbors besides N)
methyl_carbon_idx = None
for atom in mol.GetAtoms():
    if atom.GetSymbol() == "N":
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == "C":
                # Check if this carbon only has N and H as neighbors (methyl group)
                neighbor_symbols = [n.GetSymbol() for n in neighbor.GetNeighbors()]
                if neighbor_symbols.count("H") >= 3:  # Methyl group has 3 H's
                    methyl_carbon_idx = neighbor.GetIdx()
                    break

# Set up custom atom colors using drawOptions
draw_options = drawer.drawOptions()
draw_options.addAtomIndices = False
draw_options.addStereoAnnotation = False
draw_options.explicitMethyl = True  # Show methyl carbons

# Set oxygen (atomic number 8) to black for label and bonds using updateAtomPalette
draw_options.updateAtomPalette({8: (0.0, 0.0, 0.0)})

# Force all carbons to be explicitly labeled
# Also set explicit labels for oxygens
for atom in mol.GetAtoms():
    if atom.GetSymbol() == "C":
        atom.SetProp("atomLabel", "C")

# Create color dictionary for atoms
atom_colors = {}
highlight_atoms = []

for atom in mol.GetAtoms():
    idx = atom.GetIdx()

    if atom.GetSymbol() == "C":
        # All carbons blue (#1F77B4)
        atom_colors[idx] = (31 / 255, 119 / 255, 180 / 255)
        highlight_atoms.append(idx)
    elif atom.GetSymbol() == "N":
        # Nitrogen orange (#F28E2B)
        atom_colors[idx] = (242 / 255, 142 / 255, 43 / 255)
        highlight_atoms.append(idx)
    elif atom.GetSymbol() == "H":
        # Check if hydrogen is attached to the methyl carbon
        neighbors = atom.GetNeighbors()
        if len(neighbors) == 1 and neighbors[0].GetIdx() == methyl_carbon_idx:
            # This H is on the methyl group attached to N
            atom_colors[idx] = (148 / 255, 103 / 255, 189 / 255)  # Purple (#9467BD)
            highlight_atoms.append(idx)
    # Do not highlight O; label is set above so it will be drawn in black by default

# Set highlighting radius to 0 to avoid bond coloring
draw_options.fillHighlights = False
draw_options.highlightRadius = 0.3
draw_options.setHighlightColour = False

# Set options to reduce whitespace
draw_options.padding = 0

# Draw molecule with atom colors only
drawer.DrawMolecule(
    mol,
    highlightAtoms=highlight_atoms,
    highlightAtomColors=atom_colors,
    highlightBonds=[],
)
drawer.FinishDrawing()

# Save to file
drawer.WriteDrawingText("figures/structures/homarine_v2.png")
print("Image saved as homarine_v2.png")
