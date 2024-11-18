

### **Weekly Protocol (11.11.2024 – 17.11.2024)**  

**Git Branch**: `Gene-Families-tests-Andre`  

---

#### **Overview**  
This week focused on refining and expanding analysis pipelines for gene families using *Drosophila* expression data. Key tasks included developing unit tests, incorporating additional gene mapping information, and addressing gaps in species coverage in expression profiles. Significant progress was made in restructuring datasets, improving R scripts, and generating preliminary results for two dietary conditions (Diet M and Diet P).  

---

#### **Key Tasks & Progress**  

1. **Unit Testing**:  
   - Developed and documented test cases for functions in `compute_funks.R` and `load_data_funks.R` using the `testthat` package.  
   - Created a `tests` directory with initial results stored in `tests/test_results`.  
   - Paused further testing to prioritize data analysis tasks as per team discussions.  

2. **Data Preparation**:  
   - Enhanced the dataset using spreadsheets from the Diet paper, focusing on dmel and dsec gene expression data.  
   - Extracted ~2000 genes from a possible 13,000 using mapping data.  
   - Mapped dmel gene identifiers to dsec-specific identifiers to expand species coverage.  

3. **Expression Profiles**:  
   - Created four expression profiles:  
     - Two for specific diets (Diet M and Diet P).  
     - Two species-specific profiles (dmel and dsec).  
   - Combined dmel and dsec datasets into a unified format, incorporating FBgn and FBpp identifiers for compatibility.  

4. **R Script Adjustments**:  
   - Updated R scripts to handle the new data structure and address issues during result generation.  
   - Introduced log-transformed values for better analysis of metrics like Euclidean distances.  

5. **Result Generation**:  
   - Generated boxplots and CSV outputs for two test directories: `test_diet_M` and `test_diet_P`.  
   - Results now include additional categories like con_orthologs, alongside previously analyzed *special_in_paralogs*.  

6. **Discussions & Planning**:  
   - Discussed data limitations, especially the incomplete mapping for dmel to other species.  
   - Evaluated potential alternative workflows, including raw reads with kallisto.  
   - Explored online databases for additional data sources to enhance coverage.  

---

#### **Doubts & Issues**  
   - Discrepancy between significance levels in boxplots and CSV outputs—potential issue in t-test implementation.  
   - Should boxplots differentiate by species or focus solely on gene groups?  
   - Limited gene mapping success for non-dmel species—only 2000 of 13,000 genes mapped.  
   - Is the current dataset sufficient, or should alternative sources like kallisto be prioritized?  

---

#### **Next Steps**  
   - Verify and correct t-test implementation in R scripts.  
   - Expand mapping data coverage by exploring additional databases or using kallisto for raw reads.  
   - Continue refining analysis pipelines and reviewing relevant literature for alternative approaches.  
   - Investigate incorporating additional species into the analysis.  