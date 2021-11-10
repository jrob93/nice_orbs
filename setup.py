import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='nice-orbs',
    version='0.0.1',
    author='James Robinson',
    author_email='jrobinson72@qub.ac.uk',
    description='Keplerian orbit tools',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/jrob93/nice_orbs',
    # project_urls = {
    #     "Bug Tracker": "https://github.com/mike-huls/toolbox/issues"
    # },
    license='MIT',
    packages=['nice_orbs'],
    install_requires=['requests'],
)
