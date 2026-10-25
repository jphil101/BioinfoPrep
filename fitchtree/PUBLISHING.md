# Publishing `fitchtree` to PyPI

Run these on your Mac from the `fitchtree` folder of this repo:

```bash
cd "/Users/jphil/Documents/#0 - Mini-Kamikaze/#0 - Bioinfo/CodePrep/fitchtree"
```

## 0. Before the first upload

1. **Check the name is still free.** Open https://pypi.org/project/fitchtree/. It should say "404 Not Found".
2. **Check the author line** in `pyproject.toml` (`authors = [{ name = "jphil101" }]`). Everything in that
   file is shown publicly on PyPI. To add an email: `authors = [{ name = "Your Name", email = "you@example.com" }]`.
3. **Create two accounts** (they are separate sites):
   - https://test.pypi.org/account/register/ (a practice site)
   - https://pypi.org/account/register/ (the real one)

   Verify your email on both and turn on two-factor authentication. Uploads need 2FA.
4. **Create an API token on each site.** Go to Account settings → API tokens → *Add API token*, with scope
   "Entire account". Project-only scope isn't possible until the project exists. Copy the token, which starts
   with `pypi-`. It's shown only once.

## 1. Build

```bash
python3 -m pip install --upgrade build twine
rm -rf dist
python3 -m build              # creates dist/fitchtree-0.1.0.tar.gz and dist/fitchtree-0.1.0-py3-none-any.whl
python3 -m twine check dist/* # both should say PASSED
```

## 2. Practice upload to TestPyPI

```bash
python3 -m twine upload --repository testpypi dist/*
#   username: __token__
#   password: <your TestPyPI token>
```

Check that it installs and works:

```bash
python3 -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ fitchtree
python3 -c "from fitchtree import fitch_margoliash; print(fitch_margoliash([[0,22,39,39],[22,0,41,41],[39,41,0,18],[39,41,18,0]], 'ABCD').format('newick'))"
#   expected: (A:10,B:12,(C:9,D:9):20):0;
python3 -m pip uninstall -y fitchtree
```

## 3. Real upload to PyPI

```bash
python3 -m twine upload dist/*
#   username: __token__
#   password: <your PyPI token>
```

The package is live at https://pypi.org/project/fitchtree/ and `pip install fitchtree` works for everyone.
The scripts' auto-install line can then install it too.

## 4. Afterwards

- **Swap the token:** on pypi.org, delete the "Entire account" token and create one scoped to `fitchtree` only.
- **Update `requirements.txt`** at the repo root: replace the `./fitchtree` line with `fitchtree`.
- **To publish a fix later:** raise `__version__` in `src/fitchtree/__init__.py` (e.g. `0.1.1`), then repeat
  steps 1 and 3. PyPI never accepts the same version number twice, even after you delete a release.

## Optional: skip typing the token every time

Create `~/.pypirc`:

```ini
[pypi]
username = __token__
password = pypi-...your project-scoped token...
```

Then run `chmod 600 ~/.pypirc` so other users on the Mac can't read it.
