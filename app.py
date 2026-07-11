import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import tempfile
import os

st.set_page_config(page_title="AS-MS Pooling Engine", page_icon="🧪", layout="wide")

st.title("NCATS ASMS 384-Well Pooling Engine")
st.markdown("Created by Colin Gordy for the creation of pooled compound libraries for HRMS analysis. This tool supports the development of an NCATS HTS screening assay utilizing His-tagged proteins bound to Ni-NTA magnetic beads. Upload Spotfire `.sdf` library export with the NCGC ID and SMILES to predict ionization modes, maximize $m/z$ diversity, and generate 384-well plate maps for CoMa source plate prep.")

# 1. File Uploader
uploaded_file = st.file_uploader("Choose an SDF file", type=["sdf"])

def process_sdf(file_path):
    supplier = Chem.SDMolSupplier(file_path)
    data = []
    
    # Structural SMARTS patterns for ionization
    basic_nitrogen = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O);!$(N-[#6a])]")
    acidic_group = Chem.MolFromSmarts("[C,S](=[O,S])[O;H1,-1]")
    
    for idx, mol in enumerate(supplier):
        if mol is None: continue
        
        # Grab internal SAMPLE_ID
        sample_id = None
        for prop_name in ['SAMPLE_ID', 'Name', 'ID', 'sample_id', 'id']:
            if mol.HasProp(prop_name):
                sample_id = mol.GetProp(prop_name)
                break
        if not sample_id: sample_id = f"UNKNOWN_{idx}"
            
        try:
            exact_mass = Descriptors.ExactMolWt(mol)
            has_base = mol.HasSubstructMatch(basic_nitrogen)
            has_acid = mol.HasSubstructMatch(acidic_group)
            
            if has_base and not has_acid:
                predicted_mode = "positive"
                target_mz = exact_mass + 1.0073
            elif has_acid and not has_base:
                predicted_mode = "negative"
                target_mz = exact_mass - 1.0073
            else:
                predicted_mode = "positive"  # Small molecule default
                target_mz = exact_mass + 1.0073
                
            data.append({
                'SAMPLE_ID': sample_id,
                'Exact_Mass': exact_mass,
                'Target_m_z': target_mz,
                'Predicted_Mode': predicted_mode,
                'SMILES': Chem.MolToSmiles(mol)
            })
        except:
            continue
    return pd.DataFrame(data)

def assign_384_wells(df, pool_size=10):
    pooled_records = []
    
    # 384-well plate geometry mapping
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    columns = range(1, 25)
    well_coordinates = [f"{r}{c:02d}" for r in rows for c in columns]
    
    current_plate = 1
    well_pointer = 0
    
    # Process Positive and Negative modes separately to avoid polarity flipping
    for mode, group in df.groupby('Predicted_Mode'):
        sorted_group = group.sort_values(by='Target_m_z').reset_index(drop=True)
        num_compounds = len(sorted_group)
        num_pools = num_compounds // pool_size if num_compounds // pool_size > 0 else 1
        
        # Build empty structure array for the stride layout
        stride_pools = [[] for _ in range(num_pools)]
        for idx, row in sorted_group.iterrows():
            pool_idx = idx % num_pools
            if len(stride_pools[pool_idx]) < pool_size:
                stride_pools[pool_idx].append(row)
        
        # Map each generated pool to physical 384-well coordinates
        for pool in stride_pools:
            if not pool: continue
            if well_pointer >= len(well_coordinates):
                well_pointer = 0
                current_plate += 1
                
            assigned_well = well_coordinates[well_pointer]
            well_pointer += 1
            
            for sub_idx, comp in enumerate(pool):
                pooled_records.append({
                    'Destination_Plate': f"ASMS_{mode.upper()}_PLT_{current_plate}",
                    'Destination_Well': assigned_well,
                    'Well_Sub_Index': sub_idx + 1,
                    'SAMPLE_ID': comp['SAMPLE_ID'],
                    'Exact_Mass': round(comp['Exact_Mass'], 4),
                    'Target_m_z': round(comp['Target_m_z'], 4),
                    'Ionization_Mode': comp['Predicted_Mode'],
                    'SMILES': comp['SMILES']
                })
                
    return pd.DataFrame(pooled_records)

# ==========================================
# 2. Main App Execution
# ==========================================
if uploaded_file is not None:
    with st.spinner("Parsing chemical structures and generating optimal layout..."):
        
        # ✨ THE LOCAL FIX: Write the file directly inside your project folder
        local_tmp_path = "temp_upload_data.sdf"
        with open(local_tmp_path, "wb") as f:
            f.write(uploaded_file.getvalue())
            
        # Parse the local file cleanly
        raw_df = process_sdf(local_tmp_path)
        
        # Safe cleanup: Delete the local temp file immediately after reading it
        if os.path.exists(local_tmp_path):
            os.remove(local_tmp_path)
        
        if not raw_df.empty:
            # Generate the pools using the default internal names
            final_map = assign_384_wells(raw_df, pool_size=10)
            
            # Sort the data using internal names before renaming them
            final_map = final_map.sort_values(by=['Destination_Plate', 'Target_m_z']).reset_index(drop=True)
            
            # Keep your custom column names working perfectly
            final_map = final_map.rename(columns={
                'SAMPLE_ID': 'NCGC_ID',
                'Target_m_z': 'Target_M/Z'
            })
            
            # Metrics Row
            st.success(f"Successfully processed {len(raw_df)} compounds!")
            col1, col2, col3 = st.columns(3)
            col1.metric("Total Compounds", len(raw_df))
            col2.metric("Total 384-Well Pools Created", len(final_map['Destination_Well'].unique()))
            col3.metric("Total Plates Required", len(final_map['Destination_Plate'].unique()))
            
            # 3. Download Section
            st.markdown("### Download Plate Map")
            csv_data = final_map.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="Download Compound Management Manifest (CSV)",
                data=csv_data,
                file_name="asms_384well_pooling_manifest.csv",
                mime="text/csv",
            )
            
            # 4. Interactive Preview
            st.markdown("### Plate Layout Preview")
            st.dataframe(
                final_map[['Destination_Plate', 'Destination_Well', 'Well_Sub_Index', 'NCGC_ID', 'Target_M/Z', 'Ionization_Mode']], 
                use_container_width=True
            )
        else:
            st.error("Could not parse any valid structures from the uploaded SDF file.")
