# ALS 2016/12/06
# 

"""
Test that makebatch can make batch and load objects

This is not a unit test. 
"""

import unittest

import os
import shutil
import astropy.table as at
import filecmp

import bubbleimg.make_batch as make_batch

objname = 'SDSSJ1000+1242/'

# class TestCaseDownloadObj(unittest.TestCase):


#     def __init__(self, *args, **kwargs):
#         super(TestCaseDownloadObj, self).__init__(*args, **kwargs)
#         self.assertions = assertions(self)

#     @classmethod
#     def setUpClass(cls):
#         # setting up paths
#         cls.dir_test = 'bubbleimg/test/'
#         cls.dir_sampletemp = cls.dir_test+'temp_sample/'
#         cls.dir_objmock = cls.dir_test+'mock_sample/'+objname
#         cls.dir_objtemp = cls.dir_sampletemp+objname

#         # copy sample list from mock_sample/ to temp_sample/ for testing
#         if not os.path.isdir(cls.dir_sampletemp):
#             os.mkdir(cls.dir_sampletemp)

#         # copy list file to temp_sample
#         filefrom = cls.dir_test+'mock_sample/'+'list.txt'
#         fileto = cls.dir_sampletemp+'list.txt'
#         shutil.copyfile(filefrom, fileto)

#         # call make batch
#         list_torun = at.Table.read(cls.dir_sampletemp+'list.txt', format='ascii')
#         make_batch.make_batch(dir_batch=cls.dir_sampletemp, list_torun=list_torun, bandline='r', bandconti='z', catalog='mullaney', tosummarize=True)


#     @classmethod
#     def tearDownClass(cls):
#         # delete the directory 'bubbleimg/test/temp_sample/'
#         shutil.rmtree(cls.dir_sampletemp)


#     def test_create_obj_folder(self):
#         # check if obj folder is created
#         assert (os.path.isdir(self.dir_objtemp))


#     def test_download_xid(self):
#         self.assertions.assert_tempfile_equalto_mockfile(filename='xid.csv')

#     def test_download_PhotoObj(self):
#         self.assertions.assert_tempfile_equalto_mockfile(filename='PhotoObj.csv')

#     def test_download_stamp(self):
#         self.assertions.assert_tempfile_equalto_mockfile(filename='stamp-r.fits')

#     def test_download_spec(self):
#         self.assertions.assert_tempfile_equalto_mockfile(filename='spec.fits')





class TestCaseMakeLineMaps(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestCaseMakeLineMaps, self).__init__(*args, **kwargs)
        self.assertions = assertions(self)

    @classmethod
    def setUpClass(cls):
        # setting up paths
        cls.dir_test = 'bubbleimg/test/'
        cls.dir_sampletemp = cls.dir_test+'temp_sample/'
        cls.dir_objmock = cls.dir_test+'mock_sample/'+objname
        cls.dir_objtemp = cls.dir_sampletemp+objname

        # copy sample list from mock_sample/ to temp_sample/ for testing
        if not os.path.isdir(cls.dir_sampletemp):
            os.mkdir(cls.dir_sampletemp)
        if not os.path.isdir(cls.dir_objtemp):
            os.mkdir(cls.dir_objtemp)

        # copy list file to temp_sample
        filefrom = cls.dir_test+'mock_sample/'+'list.txt'
        fileto = cls.dir_sampletemp+'list.txt'
        shutil.copyfile(filefrom, fileto)
        list_torun = at.Table.read(cls.dir_sampletemp+'list.txt', format='ascii')

        # copy downloaded files
        filestocopy = ['PhotoObj.csv', 'xid.csv', 'spec.fits', 'stamp-g.fits', 'stamp-i.fits', 'stamp-r.fits', 'stamp-u.fits', 'stamp-z.fits', ]
        for filetocopy in filestocopy:
            filefrom = cls.dir_objmock+filetocopy
            fileto = cls.dir_objtemp+filetocopy
            shutil.copyfile(filefrom, fileto)

        # call make batch
        make_batch.make_batch(dir_batch=cls.dir_sampletemp, list_torun=list_torun, bandline='r', bandconti='z', catalog='mullaney', tosummarize=True)


    @classmethod
    def tearDownClass(cls):
        # delete the directory 'bubbleimg/test/temp_sample/'
        shutil.rmtree(cls.dir_sampletemp)

    def test_OIII_F_exists(self):
        filename = self.dir_objtemp+'stamp-lOIII5008_F.fits'
        assert os.path.isfile(filename)

    def test_OIII_I_exists(self):
        filename = self.dir_objtemp+'stamp-lOIII5008_I.fits'
        assert os.path.isfile(filename)


class assertions(object):

    def __init__(self, outer):
        self.outer = outer

    def assert_tempfile_equalto_mockfile(self, filename):
        # check file exist
        assert (os.path.isfile(self.outer.dir_objtemp+filename))
        # check identical to the mock one
        filepath1 = self.outer.dir_objmock+filename
        filepath2 = self.outer.dir_objtemp+filename
        assert filecmp.cmp(filepath1, filepath2)



if __name__ == '__main__':
    unittest.main()