# dtpcrm 0.1.2

* Xiaoran Lai takes over as maintainer; Christina Yap remains an author.
* `applied_titecrm_sim()` now simulates a DLT time for each patient with a
  toxicity, drawn uniformly over the observation window. A DLT only counts
  towards the model once the patient's follow-up has reached that time;
  previously all toxicities were counted as soon as the patient was dosed.
* Fixed the "Lost braces" NOTE in the `applied_crm()` documentation.
* Removed `LazyData` from DESCRIPTION as the package has no data directory.

# dtpcrm 0.1.1

* Enhanced table outputs for simulations and DTPs.
