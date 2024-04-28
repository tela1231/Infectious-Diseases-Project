import csv
from datetime import datetime

def parse_csv(csv_name):
    data = []
    with open(csv_name, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        header = next(csv_reader)  # Read the header row
        first_date = None
        for row in csv_reader:
            date_str = row[0]
            value_str = row[1]
            date = datetime.strptime(date_str, '%Y%m%d')  # Parse date string
            if first_date is None:
                first_date = date
                t = 0
            else:
                t = (date - first_date).days  # Calculate days from first date
            value = float(value_str)  # Convert value to float
            data.append((t, value))  # Append (t, value) tuple to data list
    return header, data