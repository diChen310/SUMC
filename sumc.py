#!/usr/bin/env python
# coding: utf-8


import torch
import numpy as np
import pandas as pd
from torch import nn, optim
from torch.utils.data import TensorDataset, DataLoader
from numpy.linalg import svd
from torch.nn.utils.parametrizations import orthogonal

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from tqdm import tqdm
import argparse
import os
import sys

class MLP(nn.Module):
    def __init__(self, **kwargs):
        super().__init__()
        self.layer1 = nn.Linear(kwargs["input_shape"],kwargs["hidden_shape"])
        self.layer2 = orthogonal(nn.Linear(kwargs["hidden_shape"], kwargs["hidden_shape"])) 
        self.layer3 = nn.Linear(kwargs["hidden_shape"],kwargs["input_shape"])

    def forward(self, x):
        x = self.layer1(x)
        x = self.layer2(x)
        x = self.layer3(x)
        return x

    def get_embeddings(self,x):
        x = self.layer2(self.layer1(x))
        return x


class SUMC(nn.Module):
    def __init__(self, inputs,adj_file,device='cpu',N=2000):
        super(SUMC, self).__init__()
        self.device = device
        self.inputs = torch.Tensor(inputs).to(self.device)
       
        adj = pd.read_csv(adj_file,index_col=0).to_numpy()
        
        pcs = PCA(64).fit_transform(inputs)
        self.pcs = torch.Tensor(pcs).to(self.device) 
        self.adj = torch.Tensor(adj).to(self.device)

        self.ind_views = [0,1]
        self.combinations = [(0, 1)]
        self.N = N
    def train(self):
        self.mlps = MLP(input_shape = self.pcs.shape[1], hidden_shape = 16).to(self.device) 
        def sc_loss(A,Y):
            return (torch.triu(torch.cdist(Y,Y))*torch.triu(A)).mean()
          
        #for i in range(self.num_views):
        self.mlps.train()
        optimizer = optim.Adam(self.mlps.parameters(), lr=1e-3)          
        for epoch in tqdm(range(self.N), desc='Training network '):
            optimizer.zero_grad()
            x_hat = self.mlps(self.pcs)
            Y1 = self.mlps.get_embeddings(self.pcs)
            loss1 = nn.MSELoss()(self.pcs,x_hat)
            loss2 = sc_loss(self.adj, Y1)
            loss=loss1+loss2
            loss.backward()
            optimizer.step()

        self.mlps.eval() 
        Y = self.mlps.get_embeddings(self.pcs)
        Y = StandardScaler().fit_transform(Y.cpu().detach().numpy())
        
        self.emb = Y
def main():
    """main：Process the command-line parameters and run the SUMC model"""
    parser = argparse.ArgumentParser(
        description='SUMC, Spatially Unified Meta-Clustering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python sumc.py matrix.csv adjacency.csv
  python sumc.py matrix.csv adjacency.csv -o output.csv -n 1000
  python sumc.py matrix.csv adjacency.csv --device cuda
        '''
    )
    
    # parameters
    parser.add_argument('matrix_file', help='Input the merged gene signature matrix (.csv)')
    parser.add_argument('adjacency_file', help='Input the merged similarity matrix (.csv)')
    
    # optional
    parser.add_argument('-o', '--output', default='SUMC_output_emb.csv',
                       help='Output the path of the embedded file(default: SUMC_output_emb.csv)')
    parser.add_argument('-n', '--epochs', type=int, default=1000,
                       help='Epochs (default: 1000)')
    parser.add_argument('--device', choices=['cpu', 'cuda'], default='cpu',
                       help='Device (default: cpu)')
    parser.add_argument('--hidden_shape', type=int, default=16,
                       help='Hidden shape (default: 16)')
    
    args = parser.parse_args()
    
    # check file
    if not os.path.exists(args.matrix_file):
        print(f"Error: gene signature matrix file '{args.matrix_file}' absent")
        sys.exit(1)
        
    if not os.path.exists(args.adjacency_file):
        print(f"Error: similarity matrix file '{args.adjacency_file}' absent")
        sys.exit(1)
    
    # Check if the output directory exists. If it doesn't, create it
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Create an output directory: {output_dir}")
    
    print("Read data...")
    # read data
    mat_1_df = pd.read_csv(args.matrix_file, index_col=0)
    row_names_1 = mat_1_df.index.tolist()
    mat_1 = mat_1_df.to_numpy()
    
    print(f"Read data finished:")
    print(f"  - signature shape: {mat_1.shape}")
    print(f"  - local cluster number: {len(row_names_1)}")
    print(f"  - Epochs: {args.epochs}")
    print(f"  - Device: {args.device}")
    
    # Device
    device = args.device
    if args.device == 'cuda' and not torch.cuda.is_available():
        print("Warning: CUDA is unavailable; use CPU instead")
        device = 'cpu'
    
    # Train
    print("Start training the SUMC model...")
    model = SUMC(mat_1, args.adjacency_file, device=device, N=args.epochs)
    model.train()
    
    # Obtain the embedding result
    emb_out = model.emb
    
    # Save result
    print(f"Save the embedded result to: {args.output}")
    df = pd.DataFrame(emb_out, index=row_names_1)
    df.to_csv(args.output)
    
    print("SUMC training is completed！")
    print(f"Embedded matrix shape: {emb_out.shape}")
    print(f"Output file: {args.output}")


if __name__ == "__main__":
    main()

# mat_1_df = pd.read_csv(f"/share/dichen/lungST/variables/bestk_input_matrix_0.csv",index_col=0)
# adj_file = f"/share/dichen/lungST/variables/bestk_adj_matrix_0.csv"

# row_names_1 = mat_1_df.index.tolist()
# mat_1 = mat_1_df.to_numpy()
# model = SUMC(mat_1,adj_file,N=1000)

# model.train()
# emb_out = model.emb
# df = pd.DataFrame(emb_out)
# df.to_csv('/share/dichen/lungST/pythonOut/Simcc_output_bestK_all_emb.csv')




