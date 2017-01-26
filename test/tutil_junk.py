# tutil.py

"""
store methods that will be used in test cases
"""

import filecmp

class assertions(object):

    def assert_tempfile_equalto_mockfile(self, filename):
        # check file exist
        assert (os.path.isfile(self.dir_objtemp+filename))
        # check identical to the mock one
        filepath1 = self.dir_objmock+filename
        filepath2 = self.dir_objtemp+filename
        assert filecmp.cmp(filepath1, filepath2)

