import sys
import pkg_resources
from pkg_resources import DistributionNotFound, VersionConflict

# here, if a dependency is not met, a DistributionNotFound or VersionConflict
# exception is thrown. 
with open('./requirements.txt', 'r') as f:
    dependencies = f.readlines()
    try:
        pkg_resources.require(dependencies)
    except:
        print 'Error: Some required package is missing. See the requirements file.'
        sys.exit(1)
