from argparse import ArgumentParser

#input:
#wind speed
#surface dim
#amplitude
#wave size

#class as model of input command variables evaluator
class ArgProcessor():

    def __init__(self):
        parser = ArgumentParser()
        parser.add_argument("-w", "--wind_speed", help="Wind speed, for example: 1.5", type=float, required=True)
        parser.add_argument("-d", "--grid_dim", help="Grid dimension for (x, y) grid, for example: 100", type=int, required=True)
        parser.add_argument("-a", "--amplitude", help="Default amplitude for Gerstner wave and it's potential, for example: 0.5", type=float, required=True)
        parser.add_argument("-l", "--min_wave_size", help="Minimum wave size, for example: 0.006", type=float, required=True)
        self.args = vars(parser.parse_args())

    def getInputs(self):
        return self.args
