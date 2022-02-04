"""
BASIL GUI for Oxford ASL - Base classes and functions

Copyright (C) 2020 University of Oxford
"""
import os
import functools

import nibabel as nib

class OptionError(RuntimeError):
    """
    Exception thrown because of an invalid option value (i.e.
    an error that the user must resolve)
    """

class OptionComponent():
    """
    Base class for anything which defines / responds to options
    """
    def __init__(self, app):
        self.app = app

    def _check_exists(self, label, fname, can_be_none=True):
        """
        Check if a file exists

        :param label: Human readable label for file (e.g. 'Calibration data')
        :param fname: Filename
        :param can_be_none: If None is an acceptable value

        :raises OptionError: if file does not exist
        """
        if fname is None and can_be_none:
            return
        elif not fname:
            raise OptionError("%s must be specified" % label)
        elif not os.path.exists(fname):
            raise OptionError("%s - no such file or directory" % label)

    @functools.lru_cache(maxsize=128)
    def _check_image(self, label, fname, can_be_none=True):
        """
        :param label: Human readable label for file (e.g. 'Calibration data')
        :param fname: Filename
        :param can_be_none: If None is an acceptable value

        :raises OptionError: if file does not contain an image
        """
        if fname is None and can_be_none:
            return

        self._check_exists(label, fname, can_be_none)
        try:
            nii = nib.load(fname)
        except nib.filebasedimages.ImageFileError:
            raise OptionError("%s - failed to load file - is this a valid image file?" % label)

        try:
            return nii.get_data().shape
        except:
            raise OptionError("%s - failed to read data shape - check image file is not corrupted")

    def _get_nvols(self, fname):
        """
        Get the number of volumes in a Nifti data set

        Does not throw an exception on error - use _check_image for that

        :param fname: File name
        :return Number of volumes or -1 if could not read the file for any reason
        """
        try:
            shape = self._check_image("image", fname)
            if len(shape) == 4:
                return shape[3]
            else:
                return 1
        except OptionError:
            return -1
