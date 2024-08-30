# Harmonization_methods

Propose a harmonization method to remove scanner effects(or batch effects, site effects, center effects...).

## Example of usage

Example of usage:

```
from ComBat_MLE import ComBat_MLE

#data_with_scanner_effects and scanner_labels are feature data to be harmonized and the corresponding scanner labels.

myModel=ComBat_MLE()
harmonized_data=myModel.harmonize_data(data_with_scanner_effects, scanner_labels)
```

A more detail example can be found in "test_example.ipynb".

## Features used in our experiments

We make the radiomic features used in our exeperiments available, to support future researchers to test their harmonization methods.

Radiomic features were extracted from the MRI scans of the homogeneous phantom and heterogeneous phantom, with different image preprocessing (No preprocessing, N4 bias field correction, resampling) and normalization methods (No normalization, Nyul normalization, or Z-Score normalization).  In the experiments of this MedAI 2024 paper, we take no preprocessing but only the Nyul normalization. Namely, the features used in MedAI 2024 paper are saved in the folder 'Features\No_preprocess\nyul_normalize'.

See more details about the feature extration in the [Appendix](Appendices/Appendices_to_our_MedAI_2024_paper_V1.pdf) and in our previously published paper:

> **Li Y, Ammari S, Balleyguier C, et al.**  *Impact of preprocessing and harmonization methods on the removal of scanner effects in brain MRI radiomic features[J]*. Cancers, 2021, 13(12): 3000.

## Citations

Please cite the following paper if you use the codes or the radiomic features we provided:

```
@article{li2021impact,
  title={Impact of preprocessing and harmonization methods on the removal of scanner effects in brain MRI radiomic features},
  author={Li, Yingping and Ammari, Samy and Balleyguier, Corinne and Lassau, Nathalie and Chouzenoux, Emilie},
  journal={Cancers},
  volume={13},
  number={12},
  pages={3000},
  year={2021},
  publisher={MDPI}
}
```

## Contact me

Email: yingpingleee@126.com

​            liyingping@xidian.edu.cn
