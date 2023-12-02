# Python modules


# Third party modules
import numpy as np
import functools

# Our modules



def requires_transform(func):
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if self.transform is None:
            raise ValueError("No transform set for {} object {}".format(type(self), self))
        return func(self, *args, **kwargs)

    return wrapper
    

class BaseTransform(object):
    """
    Base class for objects that require an affine tranform to make their days
    happy. In Vespa this is mainly the mrs_data_raw and mri_data classes.

    """
    def __init__(self, transform=None):
        self.transform = transform


    @requires_transform
    def to_scanner(self, *args):
        """
        Converts 3D position in BaseTransform space to the scanner reference space.
        Args x,y,z can be either 3 individual floats a numpy array_like with
        final dimension of size 3. Raises ValueError if no transform set.

        Parameters:
            x (float): The x coordinate of the point
            y (float): The y coordinate of the point
            z (float): The z coordinate of the point

        Returns:
            numpy.ndarray: The transformed 3d point in scanner coordinates

        """
        pos = normalise_positions_for_transform(*args)
        new_pos = np.einsum("ij,...j", self.transform, pos)

        return np.squeeze(np.asarray(new_pos))[..., 0:3]


    @requires_transform
    def from_scanner(self, *args):
        """
        Converts 3D position in MrScanner space to BaseTransform reference space.
        Args x,y,z can be either 3 individual floats or a numpy array_like with
        final dimension of size 3. Raises a ValueError if no transform is set.

        Parameters:
            x (float): The x coordinate of the point
            y (float): The y coordinate of the point
            z (float): The z coordinate of the point

        Returns:
            numpy.ndarray: The transformed 3d point in BaseTransform coordinates

        """
        pos = normalise_positions_for_transform(*args)
        new_pos = np.einsum("ij,...j",np.linalg.inv(self.transform), pos)

        return np.squeeze(np.asarray(new_pos))[..., 0:3]

        
    @property
    @requires_transform
    def voxel_size(self):
        """ returns, dimensions of a voxel as numpy ndarray """
        return np.linalg.norm(self.transform, axis=0)[0:3]

    @property
    @requires_transform
    def position(self):
        """ the center of the BaseTransform in scanner coordinates """
        return self.transform[:3, 3]

    @property
    @requires_transform
    def slice_vector(self):
        return self.transform[:3, 2] / np.linalg.norm(self.transform[:3, 2])

    @property
    @requires_transform
    def row_vector(self):
        return self.transform[:3, 0] / np.linalg.norm(self.transform[:3, 0])

    @property
    @requires_transform
    def col_vector(self):
        return self.transform[:3, 1] / np.linalg.norm(self.transform[:3, 1])

    @requires_transform
    def _closest_axis(self, target_axis):
        overlap = np.abs(np.dot(target_axis, self.transform[:3, :3]))
        return self.transform[:3, np.argmax(overlap)]

    @property
    @requires_transform
    def axial_vector(self):
        """
        Returns the image axis which is most closely aligned with the axial
        direction. The returned vector is guaranteed to point in the positive
        axial direction, even if the original volume vector is in the
        opposite direction.
        
        Returns:
            numpy.ndarray: The most axial image axis

        """
        # dot the three candidate vectors with (0, 0, 1)
        best_axis = self._closest_axis((0, 0, 1))
        norm_axis = best_axis / np.linalg.norm(best_axis)
        # work out if we need to reverse the direction
        return norm_axis if norm_axis[2] > 0 else -1 * norm_axis

    @property
    @requires_transform
    def coronal_vector(self):
        """
        Returns the image axis which is most closely aligned with the coronal
        direction. The returned vector is guaranteed to point in the positive
        coronal direction, even if the original volume vector is in the
        opposite direction.

        Returns:
            numpy.ndarray: The most coronal image axis

        """
        # dot the three candidate vectors with (0, 1, 0)
        best_axis = self._closest_axis((0, 1, 0))
        norm_axis = best_axis / np.linalg.norm(best_axis)
        return norm_axis if norm_axis[1] > 0 else -1 * norm_axis

    @property
    @requires_transform
    def sagittal_vector(self):
        """
        Returns the image axis which is most closely aligned with the sagittal
        direction. The returned vector is guaranteed to point in the positive
        sagittal direction, even if the original volume vector is in the
        opposite direction.

        Returns:
            numpy.ndarray: The most sagittal image axis

        """
        best_axis = self._closest_axis((1, 0, 0))
        norm_axis = best_axis / np.linalg.norm(best_axis)
        return norm_axis if norm_axis[0] > 0 else -1 * norm_axis

    @property
    @requires_transform
    def center(self):
        """ return center of image volume in scanner coords as numpy ndarray """
        return self.to_scanner((np.array(self.shape[::-1]) - 1) / 2)



def transformation_matrix(x_vector, y_vector, translation, spacing):
    """
    Creates a transformation matrix which will convert from a specified
    coordinate system to the scanner frame of reference.

    Parameters:
        x_vector (array): The unit vector along the space X axis in scanner coordinates
        y_vector (array): The unit vector along the space Y axis in scanner coordinates
        translation (array): The origin of the space in scanner coordinates
        spacing (float or array): The size of a space unit in scanner units

    Returns:
        matrix (array)

    """
    matrix = np.zeros((4, 4), dtype=np.float64)
    matrix[:3, 0] = x_vector
    matrix[:3, 1] = y_vector
    z_vector = np.cross(x_vector, y_vector)
    matrix[:3, 2] = z_vector
    matrix[:3, 3] = np.array(translation)
    matrix[3, 3] = 1.0

    # make sure that we can append to spacing
    spacing = list(spacing)
    while len(spacing) < 4:
        spacing.append(1.0)
    for i in range(4):
        for j in range(4):
            matrix[i, j] *= spacing[j]
    return matrix


def rotation_matrix(angle, axis):
    """
    Creates a 3x3 matrix which rotates `angle` radians around `axis`
    
    Parameters:
        angle (float): The angle in radians to rotate around the axis
        axis (array): The unit vector around which to rotate
        
    Returns:
        matrix (array):

    """
    c = np.cos(angle)
    s = np.sin(angle)
    matrix = np.zeros((3, 3))
    matrix[0, 0] = c + axis[0] ** 2 * (1 - c)
    matrix[0, 1] = axis[0] * axis[1] * (1 - c) - axis[2] * s
    matrix[0, 2] = axis[0] * axis[2] * (1 - c) + axis[1] * s
    matrix[1, 0] = axis[1] * axis[0] * (1 - c) + axis[2] * s
    matrix[1, 1] = c + axis[1] ** 2 * (1 - c)
    matrix[1, 2] = axis[1] * axis[2] * (1 - c) - axis[0] * s
    matrix[2, 0] = axis[2] * axis[0] * (1 - c) - axis[1] * s
    matrix[2, 1] = axis[2] * axis[1] * (1 - c) + axis[0] * s
    matrix[2, 2] = c + axis[2] ** 2 * (1 - c)
    return matrix


def normalise_positions_for_transform(*args):
    """
    Takes an input set of arguments which should represent some (x, y, z)
    coords to be transformed and makes sure they are in a numpy.ndarray,
    and adds a w dimension of magnitude 1.
    
    The two acceptable forms of input arguments are a single array_like
    with final dimension of 3, or three floating point arguments representing
    x, y and z. In the first case the returned array will have the same shape
    for all axes except the last, which will go from size 3 to 4, while in the
    second case the returned array will be of shape (4,).
    
    Parameters:
        args (array_like or 3 separate floats): The arguments to be processed
    
    Returns:
        numpy.ndarray: Points ready for transformation by a matrix

    """
    if len(args) == 3:
        positions = [*args, 1]
    elif len(args) == 1:
        positions = np.atleast_2d(args[0])
        w_array = np.expand_dims(np.ones(positions.shape[:-1]), axis=-1)
        positions = np.append(positions, w_array, axis=-1)
    else:
        raise ValueError("Unrecognised form for input args")

    return positions