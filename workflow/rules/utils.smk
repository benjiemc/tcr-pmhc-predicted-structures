def get_predicted_structure_names(path: str, template: str, seeds: list, column_name: str = 'index') -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry[column_name] for entry in csv.DictReader(fh)]

    return expand(template, name=indices, seed=seeds, sample=list(range(5)))
