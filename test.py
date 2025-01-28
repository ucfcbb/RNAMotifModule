import numpy as np
from scipy.spatial.distance import pdist, squareform

def calculate_motif_frequencies(motif_data):
    """
    Calculate expected frequencies for RNA structural motifs in 3D space.
    
    Parameters:
    motif_data: dict with keys:
        - coordinates: np.array of shape (n_motifs, 3) with x,y,z coordinates
        - types: list of motif types
        - volumes: dict mapping motif types to their volumes
    
    Returns:
    dict with expected frequencies and related statistics
    """
    
    def get_spatial_density(coords, volume):
        """Calculate local spatial density based on volume and distances"""
        # Calculate pairwise distances
        distances = pdist(coords)
        distance_matrix = squareform(distances)
        
        # Calculate average distance to nearest neighbors
        k = min(3, len(coords) - 1)  # Use up to 3 nearest neighbors
        nearest_distances = np.partition(distance_matrix, k, axis=1)[:, 1:k+1]
        avg_nearest_dist = np.mean(nearest_distances)
        
        # Estimate local density
        density = 1 / (avg_nearest_dist**3)
        return density * volume

    # Extract data
    coords = np.array(motif_data['coordinates'])
    motif_types = motif_data['types']
    volumes = motif_data['volumes']
    
    # Calculate individual frequencies
    type_counts = {}
    total_motifs = len(motif_types)
    
    for motif_type in set(motif_types):
        type_counts[motif_type] = motif_types.count(motif_type)
    
    individual_freqs = {
        motif_type: count/total_motifs 
        for motif_type, count in type_counts.items()
    }
    
    # Calculate spatial densities
    spatial_densities = {
        motif_type: get_spatial_density(
            coords[np.array(motif_types) == motif_type],
            volumes[motif_type]
        )
        for motif_type in set(motif_types)
    }
    
    # Calculate expected pairwise frequencies
    expected_frequencies = {}
    for type1 in set(motif_types):
        for type2 in set(motif_types):
            if type2 >= type1:  # Avoid duplicates
                # Basic probability multiplication
                basic_prob = individual_freqs[type1] * individual_freqs[type2]
                
                # Adjust for spatial effects
                spatial_factor = np.sqrt(
                    spatial_densities[type1] * spatial_densities[type2]
                )
                
                # Calculate final expected frequency
                expected_freq = basic_prob * spatial_factor
                
                # Store result
                expected_frequencies[(type1, type2)] = expected_freq
    
    # Calculate observed frequencies for comparison
    def count_proximal_pairs(coords, types, threshold=10.0):
        """Count pairs of motifs within threshold distance"""
        distances = pdist(coords)
        distance_matrix = squareform(distances)
        proximal_pairs = distance_matrix < threshold
        
        observed_frequencies = {}
        for type1 in set(types):
            for type2 in set(types):
                if type2 >= type1:
                    mask1 = np.array(types) == type1
                    mask2 = np.array(types) == type2
                    pair_counts = np.sum(
                        proximal_pairs[mask1][:, mask2]
                    ) / 2 if type1 != type2 else np.sum(
                        proximal_pairs[mask1][:, mask2]
                    ) / 2
                    total_possible = (
                        np.sum(mask1) * np.sum(mask2)
                    ) / 2 if type1 != type2 else (
                        np.sum(mask1) * (np.sum(mask1) - 1)
                    ) / 2
                    observed_frequencies[(type1, type2)] = (
                        pair_counts / total_possible if total_possible > 0 else 0
                    )
        return observed_frequencies

    observed_frequencies = count_proximal_pairs(coords, motif_types)
    
    # Calculate enrichment scores
    enrichment_scores = {
        pair: observed_frequencies[pair] / expected_frequencies[pair]
        if expected_frequencies[pair] > 0 else 0
        for pair in expected_frequencies.keys()
    }
    
    return {
        'individual_frequencies': individual_freqs,
        'spatial_densities': spatial_densities,
        'expected_frequencies': expected_frequencies,
        'observed_frequencies': observed_frequencies,
        'enrichment_scores': enrichment_scores
    }

# Example usage:
if __name__ == "__main__":
    # Sample data
    example_data = {
        'coordinates': [
            [0, 0, 0],
            [5, 5, 5],
            [10, 10, 10],
            [2, 2, 2],
            [7, 7, 7]
        ],
        'types': ['helix', 'loop', 'helix', 'bulge', 'loop'],
        'volumes': {
            'helix': 100,
            'loop': 50,
            'bulge': 30
        }
    }
    
    results = calculate_motif_frequencies(example_data)
    
    # Print results
    print("\nIndividual Frequencies:")
    for motif, freq in results['individual_frequencies'].items():
        print(f"{motif}: {freq:.3f}")
        
    print("\nExpected Pairwise Frequencies:")
    for pair, freq in results['expected_frequencies'].items():
        print(f"{pair}: {freq:.3f}")
        
    print("\nEnrichment Scores:")
    for pair, score in results['enrichment_scores'].items():
        print(f"{pair}: {score:.3f}")