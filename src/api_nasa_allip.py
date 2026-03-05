import requests
from astropy.table import QTable

def get_sentry_asteroids():
    # The default URL returns the summary list of all available Sentry objects
    url = "https://ssd-api.jpl.nasa.gov/sentry.api?all=1"
    
    try:
        # Get data
        print(f"Querying {url}...")
        response = requests.get(url)
        response.raise_for_status()  # Check for HTTP request errors
        
        # Parse JSON
        sentry_data = response.json()
        
        # The API returns the list of asteroids in the "data" key
        if "data" in sentry_data and sentry_data["data"]:
            data_list = sentry_data["data"]
            
            # Get the column names (keys) from the first asteroid record
            columns = list(data_list[0].keys())
            
            # Create a dictionary where keys are column names and values are lists of data
            data_dict = {
                col: [item.get(col) for item in data_list] 
                for col in columns
            }
            
            # Create the QTable
            my_qtable = QTable(data_dict)
            
            return my_qtable
            
        else:
            print("No data found in the response.")
            return None

    except requests.exceptions.RequestException as e:
        print(f"Network error: {e}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# --- Usage ---
# NASA_SENTRY_table = get_sentry_asteroids()

# if NASA_SENTRY_table:
#     print(f"\nSuccessfully created QTable with {len(NASA_SENTRY_table)} asteroids.")
    
#     # Display the first few rows
#     print("\nFirst 5 rows:")
#     print(NASA_SENTRY_table[:5])
    
#     # Show table info (columns and types)
#     print("\nTable Info:")
#     print(NASA_SENTRY_table.info)