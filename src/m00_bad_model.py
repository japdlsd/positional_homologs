import argh

import helpers
import grouping
import model_utils


def main(input_filename, output_directory, time_limit=300, threads=1):
    data = model_utils.load_task_with_metadata(input_filename)
    colors = data["colors"]
    cost = data["cost"]
    del data["colors"]
    del data["cost"]

    grouping_result = grouping.GroupingResult(list(range(1, 1 + len(colors))), 0, None, None)
    model_utils.save_results(output_directory, grouping_result, **data)


if __name__ == "__main__":
    argh.dispatch_command(main)
