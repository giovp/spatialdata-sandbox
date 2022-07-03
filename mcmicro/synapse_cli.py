import synapseclient 
import synapseutils 
import argparse
from pathlib import Path

def main(username, password):
    mcmicro_dir = Path().resolve() / "data"
    mcmicro_dir.mkdir(parents=True, exist_ok=True)

    syn = synapseclient.Synapse() 
    syn.login(username, password) 
    files = synapseutils.syncFromSynapse(syn, 'syn24849819', path=mcmicro_dir)
    
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='CLI for downloading the MCMICRO tonsil dataset')
    parser.add_argument('-u',
                        '--username',
                        action='store',
                        help='Synapse username',
                        type=str,
                        required=True)
    parser.add_argument('-p',
                        '--password',
                        action='store',
                        help='Password of the synape account',
                        type=str,
                        required=True)
    args = parser.parse_args()
    main(args.username, args.password)
    
