"""Sostituisce gli apostrofi usati al posto degli accenti (e', perche', piu') con le lettere
accentate (è, perché, più) nei file indicati. Strumento di sviluppo, usato una volta nella
rifinitura della repo pubblica; resta qui per mostrare come è stata fatta la conversione.

Non è una sostituzione cieca di «vocale + apostrofo»: nei testi dello strumento l'apostrofo
chiude anche le citazioni ('tissue type', 'scRNA'). Si convertono solo le parole di un elenco
esplicito, più «e'» quando non è il carattere 'e' fra apici.

Uso: python docs/sviluppo/accenti.py [--check] file1 [file2 ...]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

GRAVE = """unita piu identita gia puo cosi probabilita positivita numerosita modalita cioe liberta
stabilita variabilita onesta li da utilita sensibilita quantita parita generalita eleggibilita
usabilita sara circolarita verita valutabilita sparsita significativita si separabilita realta
proprieta priorita meta malignita eterogeneita esclusivita qualita capacita affidabilita
fattibilita possibilita necessita attivita""".split()
ACUTE = "perche ne poiche finche affinche benche nonche".split()
ACC = {"a": "à", "e": "è", "i": "ì", "o": "ò", "u": "ù", "A": "À", "E": "È", "I": "Ì", "O": "Ò", "U": "Ù"}
L = "A-Za-zÀ-ÿ"

_words = {w: w[:-1] + ACC[w[-1]] for w in GRAVE} | {w: w[:-1] + "é" for w in ACUTE}
_word_re = re.compile(rf"(?<![{L}])({'|'.join(sorted(_words, key=len, reverse=True))})\\?'(?![{L}])", re.IGNORECASE)
# «e'» isolato: non preceduto da un apice (sarebbe la stringa 'e'), salvo le elisioni (c'e', com'e')
_e_re = re.compile(rf"(?:(?<![{L}'])|(?<=[{L}]'))([eE])\\?'(?![{L}])")


def _fix_word(m: re.Match) -> str:
    w = m.group(1)
    out = _words[w.lower()]
    if w.isupper() and len(w) > 1:
        return out.upper()
    return out[0].upper() + out[1:] if w[0].isupper() else out


def convert(text: str) -> str:
    text = _word_re.sub(_fix_word, text)
    return _e_re.sub(lambda m: ACC[m.group(1)], text)


def main() -> None:
    args = sys.argv[1:]
    check = "--check" in args
    changed = 0
    for name in (a for a in args if a != "--check"):
        p = Path(name)
        old = p.read_text(encoding="utf-8")
        new = convert(old)
        if new != old:
            changed += 1
            print(f"{name}: {sum(a != b for a, b in zip(old.splitlines(), new.splitlines()))} righe")
            if not check:
                p.write_text(new, encoding="utf-8")
    print(f"{changed} file {'da convertire' if check else 'convertiti'}")
    if check and changed:
        sys.exit(1)


if __name__ == "__main__":
    main()
