import pickle
import json
import os

def convert_qtpy_to_json(qtpy_path, json_path):
    with open(qtpy_path, 'rb') as handle:
        data = pickle.load(handle)
    
    with open(json_path, 'w') as json_file:
        json.dump(data, json_file, indent=4)

def convert_directory(qtpy_dir, json_dir):
    os.makedirs(json_dir, exist_ok=True)
    for filename in os.listdir(qtpy_dir):
        if filename.endswith('.QTpy'):
            qtpy_path = os.path.join(qtpy_dir, filename)
            json_path = os.path.join(json_dir, filename.replace('.QTpy', '.json'))
            convert_qtpy_to_json(qtpy_path, json_path)

# Example usage:
convert_directory('./QTpy', './QTjson')
