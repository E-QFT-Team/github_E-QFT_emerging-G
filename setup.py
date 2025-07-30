#!/usr/bin/env python3
"""
Setup script for emergent-gravity package.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="emergent-gravity",
    version="1.0.0",
    author="Emergent Gravity Research Team",
    author_email="research@emergent-gravity.org",
    description="Emergent gravity from quantum field theory via holographic duality",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-username/emergent-gravity",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=2.12.0",
            "black>=21.0.0",
            "flake8>=3.9.0",
            "isort>=5.9.0",
        ],
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=0.5.0",
            "myst-parser>=0.15.0",
        ],
        "plotting": [
            "seaborn>=0.11.0",
            "plotly>=5.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "emergent-gravity=analysis.final_comprehensive_validation:main",
            "geff-scan=analysis.geff_rc_scan:main",
            "protocol-b=analysis.protocol_b_stabilization:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.md", "*.txt", "*.yml", "*.yaml"],
    },
    keywords=[
        "quantum field theory",
        "holographic duality", 
        "emergent gravity",
        "surface minimization",
        "commutator metric",
        "min-cut algorithms",
        "gravitational constant",
    ],
    project_urls={
        "Bug Reports": "https://github.com/your-username/emergent-gravity/issues",
        "Source": "https://github.com/your-username/emergent-gravity",
        "Documentation": "https://emergent-gravity.readthedocs.io/",
    },
)