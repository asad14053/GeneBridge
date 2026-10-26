import numpy as np
import faiss
from scipy.sparse import csr_matrix


def knn_graph(data, n_neighbors=5):
    """
    Compute a k-nearest-neighbor graph using FAISS.

    The query point itself is excluded.

    Important:
        At most n_samples - 1 non-self neighbors can exist.
        This function therefore safely clamps k for small minibatches.
    """

    data = np.asarray(
        data,
        dtype=np.float32,
    )

    if data.ndim != 2:
        raise ValueError(
            f"Expected 2D data, got shape {data.shape}"
        )

    n_samples, n_features = data.shape

    if n_samples == 0:
        return csr_matrix(
            (0, 0),
            dtype=np.float32,
        )

    # A one-cell minibatch has no non-self neighbors.
    if n_samples == 1:
        return csr_matrix(
            (1, 1),
            dtype=np.float32,
        )

    requested_k = int(n_neighbors)

    if requested_k <= 0:
        return csr_matrix(
            (n_samples, n_samples),
            dtype=np.float32,
        )

    # Maximum possible number of NON-SELF neighbors.
    k = min(
        requested_k,
        n_samples - 1,
    )

    index = faiss.IndexFlatL2(
        n_features
    )

    index.add(
        data
    )

    # k neighbors + the point itself.
    distances, indices = index.search(
        data,
        k + 1,
    )

    # Remove self-neighbor.
    neighbor_indices = (
        indices[:, 1:k + 1]
    )

    neighbor_distances = np.sqrt(
        distances[:, 1:k + 1]
    )

    # Defensive protection against invalid FAISS indices.
    valid = (
        neighbor_indices >= 0
    )

    rows = np.repeat(
        np.arange(
            n_samples
        ),
        k,
    )

    cols = (
        neighbor_indices
        .reshape(-1)
    )

    values = (
        neighbor_distances
        .reshape(-1)
        .astype(
            np.float32,
            copy=False,
        )
    )

    valid_flat = (
        valid.reshape(-1)
    )

    rows = rows[
        valid_flat
    ]

    cols = cols[
        valid_flat
    ]

    values = values[
        valid_flat
    ]

    sparse_matrix = csr_matrix(
        (
            values,
            (
                rows,
                cols,
            ),
        ),
        shape=(
            n_samples,
            n_samples,
        ),
        dtype=np.float32,
    )

    return sparse_matrix
