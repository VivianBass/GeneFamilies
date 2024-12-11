
**Date**: 14.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`  

---

### **Tasks**  

1. **Test Directories**:  
   - Created two new test directories: `test_diet_M` and `test_diet_P`.  
   - Focused on *Drosophila melanogaster* (dmel) and *Drosophila sechellia* (dsec) using additional data from the spreadsheets linked in the Diet paper.  
   - Addressed the issue where expression profiles contained only dmel gene identifiers, requiring mapping to dsec-specific identifiers to enable results for both species.  

2. **Gene Mapping and Expression Profiles**:  
   - Extracted ~2000 genes from 13,000 in the Diet paper mapping table that overlapped with the expression profiles.  
   - Filtered the expression profiles (covering 3 diets, 4 tissues, and 4 species) for dmel and dsec data.  
   - Created four expression profiles:  
     - Two diet-specific profiles (Diet M and Diet P).  
     - Separate dmel and dsec profiles with newly mapped dsec gene identifiers.  
   - Combined dmel and dsec profiles into a single dataset by appending dsec data below dmel data.  
   - Addressed an issue where gene groups used FBpp (protein sequence) identifiers, while expression profiles used FBgn (nucleotide sequence) identifiers:  
     - Extracted mapping information from combined FASTA files for all *Drosophila* species, linking FBgn to FBpp identifiers.  
     - Updated expression profiles to include both FBgn and FBpp identifiers for compatibility with gene groups.  
   - Finalized two expression profiles (Diet M and Diet P) for both dmel and dsec, including data for four tissues: whole body, muscle, fat, and gut.  
1. **Gene Identifier Discussion**:  
   - Discussed whether to use FBpp or FBgn identifiers.  
   - Decided to include both in the analysis by leveraging a mapping table generated from the FASTA file.  
   - Extracted parent FBgn and the longest-sequence FBpp identifiers from the FASTA file.  
---

### **Doubts and Issues**  

- Are we correctly mapping gene identifiers (FBgn to FBpp) for consistency with gene groups?  
- Can current data reliably support mapping across all species?  

---

### **Next Steps**  

