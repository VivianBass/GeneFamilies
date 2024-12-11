

### **Weekly Protocol (04.11.2024 – 08.11.2024)**

**Git Branch**: `Gene-Families-tests-Andre`

---

#### **Overview**
This week focused on refining R scripts to align with new data formats for gene groups and expression profiles, troubleshooting data intersection issues, and enhancing statistical and plotting functions for gene family experiments.

#### **Key Tasks & Progress**
- **Script Adjustments**: Updated all R scripts to handle new data formats, switching from `FBgn` gene IDs to `FBpp` protein sequences for expression profiles. Verified these updates across multiple test directories.
- **Expression Profile Testing**: Tested with alternative expression profiles from the diet paper and FlyBase to find sufficient intersection with gene group data, though results showed limited overlap, especially for `con-orthologs` and `out_paralogs`.
- **Data Filtering**: Refined dataframes to retain only intersecting genes to reduce NA values, improve computation speed by reducing unnecessary congestion and generate meaningful plots. Created two filtering approaches (v1 and v2) and confirmed they yield consistent results.
- **Distance & Statistical Calculations**: Enhanced Euclidean distance calculations, especially by adding logarithmic distance calculations. Improved significance testing for plots using `geom_signif()` from `ggsignif` and created a function to denote significance levels (`***`, `**`, `*`, or `ns`).
- **Plotting**: Streamlined plotting scripts for generating visuals in each results folder.

#### **Doubts & Issues**
   - **Data Intersection**: Limited data overlap with expression profiles, especially for `con-orthologs`, makes certain analyses challenging. Only `specific in-paralogs` and `in-paralogs` have sufficient data.
   - **Euclidean Distance for Families**: Uncertainty about including gene families in Euclidean distance calculations.
   - **T-test Accuracy**: T-test results do not perfectly align with significance levels in plots, as `geom_signif()` directly incorporates significance without using t-test results.

#### **Next Steps**
   - Integrate logarithmic distance calculations into scripts.
   - Finalize unit tests for all functions with varied test scenarios.
   - Adjust t-test result calculations to ensure alignment with plot significance levels.