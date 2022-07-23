import synapseclient 
import synapseutils 
import argparse
from pathlib import Path

def main(username, password):
    wsi_codex_dir = Path().resolve() / "data" / "WSI_tonsil" / "CODEX"
    wsi_codex_dir.mkdir(parents=True, exist_ok=True)

    syn = synapseclient.Synapse() 
    syn.login(username, password) 
    files = synapseutils.syncFromSynapse(syn, 'syn25173958', path=wsi_codex_dir)
    # files = synapseutils.syncFromSynapse(syn, 'syn22345748', path=tma_dir)


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
    
