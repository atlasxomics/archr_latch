from slims.slims import Slims

from wf.creds import username, password

def get_pk(cntn_id, slims):
    return slims.fetch('Content', equals('cntn_id', cntn_id))[0].pk()  

def push_result(payload, slims):
    return slims.add("Result", payload)

def slims_init(username=username, password=password):
    return Slims("slims", "https://slims.atlasxomics.com/slimsrest", username, password) 
