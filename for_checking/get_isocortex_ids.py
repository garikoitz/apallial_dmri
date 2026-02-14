#!/usr/bin/env python3
"""Get all Allen CCF structure IDs for ROIs needed for tractography.

ROIs needed:
(1) Isocortex
(2) Piriform cortex + Entorhinal area + Piriform amygdalar area + Postpiriform transition area
(3) Striatum dorsal region + Striatum ventral region
(4) Thalamus
(5) Cerebellum
"""

import json

with open('local/atlas/allen_ccf/structures.json') as f:
    data = json.load(f)

structures = data.get('msg', [])

def find_structure(structs, target_id):
    """Find a structure by ID in nested hierarchy"""
    for s in structs:
        if s['id'] == target_id:
            return s
        if 'children' in s:
            result = find_structure(s['children'], target_id)
            if result:
                return result
    return None

def find_structure_by_name(structs, name):
    """Find a structure by name (case-insensitive partial match)"""
    for s in structs:
        if name.lower() in s['name'].lower() or name.lower() in s.get('acronym', '').lower():
            return s
        if 'children' in s:
            result = find_structure_by_name(s['children'], name)
            if result:
                return result
    return None

def get_all_ids(struct):
    """Get all IDs from a structure and its descendants"""
    ids = [struct['id']]
    if 'children' in struct:
        for child in struct['children']:
            ids.extend(get_all_ids(child))
    return ids

# Define the parent structure IDs we need
# Based on Allen CCF ontology:
ROI_DEFINITIONS = {
    'isocortex': {
        'parent_id': 315,
        'description': 'Isocortex (all cortical areas)'
    },
    'piriform': {
        'parent_id': 961,
        'description': 'Piriform area'
    },
    'entorhinal': {
        'parent_id': 909,
        'description': 'Entorhinal area'
    },
    'piriform_amygdalar': {
        'parent_id': 788,
        'description': 'Piriform-amygdalar area'
    },
    'postpiriform': {
        'parent_id': 566,
        'description': 'Postpiriform transition area'
    },
    'striatum_dorsal': {
        'parent_id': 485,
        'description': 'Striatum dorsal region (Caudoputamen)'
    },
    'striatum_ventral': {
        'parent_id': 493,
        'description': 'Striatum ventral region (includes NAc)'
    },
    'thalamus': {
        'parent_id': 549,
        'description': 'Thalamus'
    },
    'cerebellum': {
        'parent_id': 512,
        'description': 'Cerebellum'
    }
}

print("=" * 80)
print("Allen CCF Structure IDs for Tractography ROIs")
print("=" * 80)

all_roi_ids = {}

for roi_name, roi_def in ROI_DEFINITIONS.items():
    parent_id = roi_def['parent_id']
    struct = find_structure(structures, parent_id)
    
    if struct:
        all_ids = get_all_ids(struct)
        all_roi_ids[roi_name] = all_ids
        print(f"\n{roi_def['description']}:")
        print(f"  Parent ID: {parent_id}")
        print(f"  Name: {struct['name']}")
        print(f"  Acronym: {struct.get('acronym', 'N/A')}")
        print(f"  Total IDs (including children): {len(all_ids)}")
        print(f"  IDs: {all_ids}")
    else:
        print(f"\n{roi_def['description']}: NOT FOUND (ID {parent_id})")

# Summary for config.yaml
print("\n" + "=" * 80)
print("SUMMARY FOR config.yaml (use parent IDs, atlas_manager gets all children)")
print("=" * 80)
print("""
rois:
  # (1) Isocortex (all cortical areas)
  isocortex:
    name: "Isocortex"
    allen_ids: [315]  # Parent ID - atlas_manager gets all 374 child IDs
  
  # (2) Piriform complex (multiple structures combined)
  piriform_complex:
    name: "Piriform_Complex"
    allen_ids: [961, 909, 788, 566]  # PIR, ENT, PAA, TR
  
  # (3) Striatum (dorsal + ventral combined)
  striatum:
    name: "Striatum"
    allen_ids: [485, 493]  # Caudoputamen (dorsal) + Striatum ventral region
  
  # (4) Thalamus
  thalamus:
    name: "Thalamus"
    allen_ids: [549]  # Thalamus parent
  
  # (5) Cerebellum
  cerebellum:
    name: "Cerebellum"
    allen_ids: [512]  # Cerebellum parent
""")
