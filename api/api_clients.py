import requests
import time
from core.models import GeneData
from core.utils import validate_strand

class EnsemblClient:
    def __init__(self):
        self.base_url = "https://rest.ensembl.org"
        self.mygene = MyGeneClient()

    def _handle_request(self, endpoint, headers):
        try:
            response = requests.get(f"{self.base_url}{endpoint}", headers=headers, timeout=15)
            if response.status_code == 200:
                return response
            elif response.status_code == 429:
                retry_after = int(response.headers.get("Retry-After", 1))
                time.sleep(retry_after)
                return self._handle_request(endpoint, headers)
            return None
        except:
            return None

    def get_sequence_by_id(self, ensembl_id, chunk_size=None):
        lookup_endpoint = f"/lookup/id/{ensembl_id}?"       
        gene_info_res = self._handle_request(lookup_endpoint, {"Content-Type": "application/json"})
        if not gene_info_res:
            return None       # metadata

        gene_info = gene_info_res.json()
        strand = gene_info.get('strand', 1)    #extract

        headers = {"Content-Type": "text/plain"}
        res = self._handle_request(f"/sequence/id/{ensembl_id}", headers)    
        
        if not res:
            return None
            
        sequence = res.text
        if chunk_size:
            sequence = sequence[:chunk_size]
            
        return {
            "sequence": sequence,    # object so we would know if we need to reverse index logic
            "strand": strand,
            "start": gene_info.get('start'),
            "end": gene_info.get('end')
        }

    def get_gene_info(self, ensembl_id):
        headers = {"Content-Type": "application/json"}
        res = self._handle_request(f"/lookup/id/{ensembl_id}", headers)
        return res.json() if res else None
    
    def get_gene_data(self, symbol: str) -> GeneData:
        ensembl_id = self.mygene.get_ensembl_id(symbol)
        if not ensembl_id:
            raise ValueError(f"Could not find Ensembl ID for symbol: {symbol}")
        
        gene_info = self.get_gene_info(ensembl_id)
        if not gene_info:
            raise ValueError(f"Could not retrieve gene info for Ensembl ID: {ensembl_id}")
        
        sequence_data = self.get_sequence_by_id(ensembl_id)
        if not sequence_data:
            raise ValueError(f"Could not retrieve sequence for Ensembl ID: {ensembl_id}")
        
        return GeneData(
            symbol=symbol,
            ensembl_id=ensembl_id,
            sequence=sequence_data['sequence'],
            strand=sequence_data['strand'],
            start=sequence_data['start'],
            end=sequence_data['end']
        )   
    
class MyGeneClient:
    def __init__(self):
        self.base_url = "https://mygene.info/v3"

    def get_ensembl_id(self, gene_symbol):
        params = {"q": f"symbol:{gene_symbol}", "species": "human", "fields": "ensembl.gene", "size": 1}
        try:
            response = requests.get(f"{self.base_url}/query", params=params, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if data.get('hits'):
                    ensembl_data = data['hits'][0].get('ensembl')
                    if isinstance(ensembl_data, list):
                        return ensembl_data[0]['gene']
                    return ensembl_data.get('gene')
        except:
            pass
        return None