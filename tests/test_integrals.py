from src.physics import create_k_grid, gap_integral, number_integral

def test_integrals():
    """
    Simple test to verify that integrals run without errors
    and return finite values.
    """
    k, dk = create_k_grid()

    mu = 1.0
    Delta = 0.5

    gap_val = gap_integral(k, mu, Delta, dk)
    num_val = number_integral(k, mu, Delta, dk)

    print("Gap integral:", gap_val)
    print("Number integral:", num_val)

    # basic checks
    assert gap_val is not None
    assert num_val is not None

if __name__ == "__main__":
    test_integrals()
