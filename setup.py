from setuptools import setup, find_packages

setup(
    name="Manual/windtunnel",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        # Core dependencies with flexible versions to avoid conflicts
        "pandas",
        "numpy", 
        "matplotlib",
        "scikit-learn",
        # Optional dependencies that might conflict with anaconda
        "numpy-stl; platform_system!='Windows'",  # Skip on Windows to avoid conflicts
    ],
    python_requires=">=3.7",
    # Jupyter notebook optimizations
    zip_safe=False,  # Prevents issues with Jupyter's import system
    include_package_data=True,  # Ensures all files are included
    # Handle conflicts gracefully
    extras_require={
        "full": [
            "numpy-stl>=2.16.0",
            "PyQt5>=5.15.0",
        ],
        "minimal": [
            # Only essential packages
        ],
    },
)
