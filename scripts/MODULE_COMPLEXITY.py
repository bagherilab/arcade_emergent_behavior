from .utilities import load, load_tar
import tarfile

class MODULE_COMPLEXITY():
    NAME = "MODULE_COMPLEXITY"

    METABOLISM = ['R','S','M','C']

    SIGNALING = ['R','S','M','C']

    @staticmethod
    def run(output_path, input_path, func, name=NAME,
            metabolism=METABOLISM, signaling=SIGNALING, timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for meta in metabolism:
            for sig in signaling:
                code = f"_{meta}_{sig}"
                infile = f"{input_path}{name}/{name}_{meta}_{sig}.pkl"
                print(f"{name} : {code}")

                loaded = load(infile)
                func(*loaded, outfile, code, exclude=[-1], timepoints=timepoints, seeds=seeds)

    @staticmethod
    def loop(output_path, func1, func2, extension, name=NAME,
             metabolism=METABOLISM, signaling=SIGNALING, timepoints=[]):
        outfile = f"{output_path}{name}/{name}"
        out = { "data": [] }
        tar = load_tar(outfile, extension)

        for t in timepoints:
            for meta in metabolism:
                for sig in signaling:
                    code = f"_{meta}_{sig}"
                    func1(outfile, out, { "time": t, "meta": meta, "sig": sig }, extension, code, tar=tar)

        func2(outfile, extension, out)

    @staticmethod
    def load(output_path, input_path, func, extension="", name=NAME,
             metabolism=METABOLISM, signaling=SIGNALING, timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for meta in metabolism:
            for sig in signaling:
                code = f"_{meta}_{sig}"
                infile = f"{input_path}{name}{extension}/{name}_{meta}_{sig}{extension}.tar.xz"
                print(f"{name} : {code}")

                tar = tarfile.open(infile)
                func(tar, timepoints, { "meta": meta, "sig": sig }, outfile, code)
