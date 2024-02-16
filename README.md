# PyRASOrbitalKit

Python package for running RASCI calculations using PyQchem. The main feature allows performing various types of calculations using natural orbitals computed from an initial RASCI calculation. This is achieved by employing different functions (defined in ClcFs.py file) that communicate the calculation type to QChem via PyQchem. The typical result from a calculation includes a QChem output and an electronic structure dictionary, both managed within PyQchem tools. The current version is designed to work exclusively within the first RASCI root; further developments for use within higher energy roots are necessary.

See examples for current calculations.

- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Instalation

This package can be install using pip:

```bash
pip install pyrasorbitalkit
```

## Usage
import pyrasorbitalkit

# Use a given calculation function
pyrasorbitalkit.ClcFs.RHF_calc(mol, bas)

# Parse information from RASNOF0 calculations
pyrasorbitalkit.parser_rasci_rasnof0.parser_rasci(output)


## Contributing

Additional functions or schemes working with other orbital sets are welcome.

Pull requests are welcome. For major changes, please open an issue first to discuss what would like to change.

For additional information, usage details or usage in Atlas Cluster, write to my e-mail.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
