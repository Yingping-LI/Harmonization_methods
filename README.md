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


## Contact me

Email: yingpingleee@126.com
