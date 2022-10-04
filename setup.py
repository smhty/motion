import setuptools
import dorna_motion
with open("README.md", "r") as fh:
    readme = fh.read()

setuptools.setup(
    name="dorna_motion",
    author="Dorna Robotics",
    version=dorna_motion.__version__,
    author_email="info@dorna.ai",
    description="Dorna motion API",
    long_description=readme,
    long_description_content_type='text/markdown',
    url="https://dorna.ai/",    
    packages=setuptools.find_packages(),
    classifiers=[
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3.10',
        "Operating System :: OS Independent",
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires=[],
    license="MIT",
    include_package_data=True,
    zip_safe = False,
)