"""PANTS offline pipeline: recall, negatives, embed, train, structure, harmonise.

Nothing in this package runs in a web request. All heavy compute is precomputed on the
M1 Max and shipped to the droplet as SQLite rows plus static mmCIF (spec section 3.1).
"""
