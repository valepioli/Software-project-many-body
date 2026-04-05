print("BEC–BCS project started")

from physics import create_k_grid, gap_integral, number_integral

def main():
    # Create momentum grid
    k, dk = create_k_grid()

    # Example initial guess
    mu = 1.0
    Delta = 0.5

    # Compute integrals
    gap_val = gap_integral(k, mu, Delta, dk)
    num_val = number_integral(k, mu, Delta, dk)

    # Print results for verification
    print("Gap integral value:", gap_val)
    print("Number integral value:", num_val)

if __name__ == "__main__":
    main()
