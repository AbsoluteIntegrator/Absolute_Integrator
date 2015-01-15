import skimage
import numpy as np


def subtract_background(image):
    """Subtracts background from image, isolating peaks from one another more clearly."""
    seed = np.copy(image)
    seed[1:-1, 1:-1] = image.min()
    mask = image
    dilated = skimage.morphology.reconstruction(seed, mask, method='dilation')
    return image-dilated