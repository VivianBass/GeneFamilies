
**Date**: 30.10.2024  
**Git Branch**: `Gene-Families-tests-Andre`

**Tasks**:

- load expression data , try if it works

- adjust minor issues with the loading Scripts

- unit tests for the functions and do a functions file for the loading script functions

- exact definition of the gene groups

- VPN, sophos

- load expression data , try if it works

- get information from vv on how the lists should look like (nested etc., 1st layer second etc.)

- see how to change dists rscript for logarythmic data


- Add all functions that we need to a special function file called "load_data_funks.R" inside R directory and delete them from your scripts. Call the file in your scripts so we don't have functions inside the exec scripts anymore.
- All the functions inside your load_data_funks.R should have unit tests. Please read:
	https://smartbear.com/learn/automated-testing/what-is-unit-testing/
	https://www.geeksforgeeks.org/unit-testing-in-r-programming/


- Delete numbers from scripts, we have the rscript_execution.md to know the order
- rscript_execution.md is not updated
- input_files.md is not updated


- chanmge Gene-groups lists to that format


about the format for the lists, maybe we can do something like this:
 
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

what do you think?


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

**Doubts and Issues**:

**Next Steps**:

---

**Code**