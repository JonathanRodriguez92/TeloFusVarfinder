import requests
from bs4 import BeautifulSoup
import pandas as pd
def extractionLongevity(url):

    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, "html.parser")
        lifespan_element = soup.find("h3", id="lifespan_longevity").find_next("dd")
        lifespan_text = lifespan_element.get_text(strip=True)
        words = lifespan_text.split()
        for word in words:
            try:
                years = float(word.strip())
                print(years)
                return years
            except ValueError:
                pass
    else:
        print("No se pudo acceder a la página web.") 
        
excel_file = "zoonomia_sp_info.xlsx"
df = pd.read_excel(excel_file)
df = 
species = df[['Species']]
for specie in species:
    url = f"https://animaldiversity.org/accounts/{specie}"
    longevity = extractionLongevity(url)
    print("Longevidad:", longevity)
