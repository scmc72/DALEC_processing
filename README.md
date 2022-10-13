Here's my lovely DALEC processing code.

Right now I'm going to use this to remind me of some useful things in case I forget.



### How to load functions:
Make sure that the functions directory is added to the path before attempting to import:

```
import sys
# append functions directory to the path: 
sys.path.append('C:/Users/daa5/Project/DALEC_processing/functions/')
# then can load a library:
import dalecLoad as dl
```

OR alternatively, for Jupyter notebooks (more portable):
```
import os
import sys
lib_path = os.path.abspath(os.path.join(os.path.abspath(''), '../../functions/'))
sys.path.append(lib_path)
# need to append our functions dir to the path! 
import dalecLoad as dl
```

AND, for .py scripts:

```
import os
import sys
lib_path = os.path.abspath(os.path.join(os.path.abspath(''), '../functions/'))
sys.path.append(lib_path)
# need to append our functions dir to the path! 
import dalecLoad as dl
```