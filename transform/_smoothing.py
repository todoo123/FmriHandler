import multiprocessing as mp
import numpy as np
from scipy.ndimage import convolve1d

def fwhm_to_sigma(fwhm_samples: float) -> float:
    """Convert Full Width at Half Maximum (FWHM) to standard deviation (sigma)."""
    return fwhm_samples / 2.354820045

def gaussian_kernel_1d(sigma_samples: float, truncate: float = 3.0) -> np.ndarray:
    """Generate a 1D Gaussian kernel."""
    if sigma_samples <= 0:
        return np.array([1.0])
    r = int(np.ceil(truncate * sigma_samples))
    x = np.arange(-r, r + 1)
    k = np.exp(-0.5 * (x / sigma_samples)**2)
    return k / np.sum(k)

def _smooth_chunk(args):
    """Helper function for chunk processing."""
    data_chunk, kernel = args
    return convolve1d(data_chunk, kernel, axis=-1, mode='reflect')

def gaussian_smoothing_1d(
    data: np.ndarray, 
    fwhm_samples: float = 3.0, 
    truncate: float = 3.0, 
    chunk_size: int = 20000,
    n_jobs: int = 1
) -> np.ndarray:
    """
    Apply 1D Gaussian kernel smoothing along the time axis (last axis).
    Handles large datasets by chunking and optional parallel processing.
    
    Args:
        data: Input fMRI data, typically shape (S, V, T) or (V, T).
        fwhm_samples: Full Width at Half Maximum in samples.
        truncate: Truncate the Gaussian kernel at this many standard deviations.
        chunk_size: Number of voxels/series to process at once to save memory.
        n_jobs: Number of parallel jobs. Set to -1 to use all CPU cores.
        
    Returns:
        np.ndarray: Smoothed data with the same shape as input.
    """
    original_shape = data.shape
    
    # Flatten to 2D (N, T) where N is all dimensions except the last (time)
    if data.ndim >= 2:
        T = data.shape[-1]
        data_2d = data.reshape(-1, T)
    else:
        raise ValueError("Data must be at least 2D (Voxels x Time)")
        
    N, T = data_2d.shape
    
    sigma = fwhm_to_sigma(fwhm_samples)
    kernel = gaussian_kernel_1d(sigma, truncate=truncate)
    
    # Pre-allocate output
    out_2d = np.empty_like(data_2d, dtype=np.float64)
    
    # Prepare chunks
    idx_starts = list(range(0, N, chunk_size))
    chunks = [(data_2d[i:i+chunk_size], kernel) for i in idx_starts]
    
    if n_jobs == -1 or n_jobs is None:
        n_jobs = mp.cpu_count()
        
    if n_jobs == 1:
        # Sequential processing
        for i, chunk_args in zip(idx_starts, chunks):
            out_2d[i:i+chunk_size] = _smooth_chunk(chunk_args)
    else:
        # Parallel processing
        with mp.Pool(processes=n_jobs) as pool:
            results = pool.map(_smooth_chunk, chunks)
            
        for i, res in zip(idx_starts, results):
            out_2d[i:i+chunk_size] = res
            
    return out_2d.reshape(original_shape)
