
**Date**: 30.10.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

- Resolved issues in `load_expression_data.R` rscript.

- Fixed minor issues in loading scripts (1. Loading Section).

- Clarified with Vivian the required format for gene group lists, 
  which should include nested information beyond gene family (Orthogroup) members.

- Updated gene group lists to the new nested format as agreed.

```R
    {
    'OG0000000': {
        ( 'dana','FBpp0117097'): {
            'dmoj': ['FBpp0172663'],
            'dpse': ['GA19158_PA'],
            'dsim': ['FBpp0313116'],
            'dvir': ['FBpp0236954'],
            'dwil': ['FBpp0241545'],
            'dyak': ['FBpp0256338']
        },
        ( 'dana','FBpp0118912'): {
            'dere': ['FBpp0131563'],
            'dmoj': ['FBpp0164096'],
            'dpse': ['GA19158_PA'],
            'dsec': ['FBpp0201147'],
            'dvir': ['FBpp0229725'],
            'dwil': ['FBpp0254305'],
            'dyak': ['FBpp0262422']
        },
        ...
    },
    ...
    }

    {
    'family_id': {
        ( 'Gene_species_1','gene'): {
            'Ortholog_species_1': ['ortholog'],
            'Ortholog_species_2': ['ortholog'],
            'Ortholog_species_3': ['ortholog'],
            'Ortholog_species_4': ['ortholog'],
            'Ortholog_species_5': ['ortholog']
        },
        ( 'Gene_species_2','gene'): {
            'Ortholog_species_1': ['ortholog'],
            'Ortholog_species_2': ['ortholog'],
            'Ortholog_species_3': ['ortholog'],
            'Ortholog_species_4': ['ortholog'],
            'Ortholog_species_5': ['ortholog']
        },
        ...
    },
    ...
    }
```

**Doubts and Issues**:

**Next Steps**:

- Create unit tests for functions if any are used right now in the Loading Section 
  and begin adjusting R scripts for statistical analyses and distance calculations.

---

**Code**