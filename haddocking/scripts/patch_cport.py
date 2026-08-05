"""
Patches applied to a freshly cloned haddocking/cport checkout, before `pip install .` runs it
into its own conda env. Needed because:

- scriber's server (biomine.cs.vcu.edu) sits behind a permanent anti-bot access wall, so it is
  dropped from the predictor list used by ProtCPORTPrediction (see predictStep). Both of cport's
  trained ML consensus classifiers require scriber, so neither can run without it; pwchem's
  protocol falls back to averaging the raw per-predictor scores instead (see parseMLOutput).
- cport crashes the whole run with a bare TypeError if any single predictor times out or fails
  (e.g. scannet's real backend is a queued service and often exceeds cport's 18 min retry
  budget) instead of just skipping that predictor's contribution.

Run as: python3 patch_cport.py <path-to-cport-checkout-root>
"""
import pathlib
import sys


def patch(path, old, new):
    text = path.read_text()
    if new in text:
        return  # already patched (e.g. a resumed/re-run install)
    if old not in text:
        raise RuntimeError(f"expected pattern not found in {path}, cport's source may have changed")
    path.write_text(text.replace(old, new, 1))


def main(cportRoot):
    src = pathlib.Path(cportRoot) / 'src' / 'cport'

    patch(src / 'url.py',
          'SCRIBER_URL = "http://biomine.cs.vcu.edu/servers/SCRIBER/"',
          'SCRIBER_URL = "https://biomine.cs.vcu.edu/servers/SCRIBER/"')
    patch(src / 'url.py',
          'SCANNET_URL = "http://bioinfo3d.cs.tau.ac.il/ScanNet/index_real.html"',
          'SCANNET_URL = "https://bioinfo3d.cs.tau.ac.il/ScanNet/index_real.html"')

    patch(src / 'modules' / 'scannet.py',
          '        browser.follow_link(browser.links()[7])\n'
          '        processing_url = browser.get_url()',
          '        result_link = browser.links()[7]\n'
          '        result_link["href"] = result_link["href"].replace("http://", "https://", 1)\n'
          '        browser.follow_link(result_link)\n'
          '        processing_url = browser.get_url()')

    utilsPy = src / 'modules' / 'utils.py'
    patch(utilsPy,
          '    standardized_dic = standardize_residues(result_dic, chain_id, pdb_file)\n'
          '    reslist = get_residue_range(standardized_dic)\n'
          '    data = []\n'
          '    for pred in result_dic:\n'
          '        row = [pred]\n'
          '\n'
          '        if pred in scored_predictors:',
          '    standardized_dic = standardize_residues(result_dic, chain_id, pdb_file)\n'
          '    reslist = get_residue_range(standardized_dic)\n'
          '    data = []\n'
          '    for pred in result_dic:\n'
          '        if result_dic[pred] is None:\n'
          '            log.warning(f"Predictor {pred} did not return a result, skipping it")\n'
          '            continue\n'
          '\n'
          '        row = [pred]\n'
          '\n'
          '        if pred in scored_predictors:')
    patch(utilsPy,
          '    active_reslist = []\n'
          '    passive_reslist = []\n'
          '    for pred in result_dic:\n'
          '        if pred in scored_predictors:',
          '    active_reslist = []\n'
          '    passive_reslist = []\n'
          '    for pred in result_dic:\n'
          '        if result_dic[pred] is None:\n'
          '            continue\n'
          '        if pred in scored_predictors:')
    patch(utilsPy,
          '    if bias != 0:\n'
          '        for pred in result_dic:\n'
          '            if pred not in scored_predictors and pred not in pdb_predictors:',
          '    if bias != 0:\n'
          '        for pred in result_dic:\n'
          '            if result_dic[pred] is None:\n'
          '                continue\n'
          '            if pred not in scored_predictors and pred not in pdb_predictors:')
    patch(utilsPy,
          '        for pred in result_dic:\n'
          '            # these use sequences and do not see missing items\n'
          '            if pred == "predictprotein" or pred == "scriber":',
          '        for pred in result_dic:\n'
          '            if result_dic[pred] is None:\n'
          '                continue\n'
          '            # these use sequences and do not see missing items\n'
          '            if pred == "predictprotein" or pred == "scriber":')


if __name__ == '__main__':
    main(sys.argv[1])
