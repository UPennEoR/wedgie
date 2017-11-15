from setuptools import setup
import glob
import os

def package_files(package_dir, subdirectory):
    # walk the input package_dir/subdirectory
    # return a package_data list
    paths = []
    directory = os.path.join(package_dir, subdirectory)
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            path = path.replace(package_dir + '/', '')
            paths.append(os.path.join(path, filename))
    return paths
data_files = package_files('wedgie', 'calibrations')

setup_args = {
    'name': 'wedgie',
    'author': 'UPenn EoR',
    'url': 'https://github.com/UPennEoR/wedgie',
    'license': 'BSD',
    'description': 'methods for generating polarization visibility wedges',
    'package_dir': {'wedgie': 'wedgie'},
    'packages': ['wedgie'],
    'include_package_data': True,
    'scripts': glob.glob('scripts/*'),
    'install_requires': ['aipy>=2.1.6', 'pyuvdata>=1.2', 'numpy', 'astropy>=1.2',
                         'matplotlib', 'hera_cal>=1.0'],
    'package_data': {'wedgie': data_files},
    'version': 0.1,
}

if __name__ == '__main__':
    apply(setup, (), setup_args)
