Here we start to use a set of sims from Padma in Mar 2026. 

We fix the RA/DEC coordinate of each stamp. 

Currently we have ~3000 lens systems. We want to include them all (or as many as possible) in one template/visit. 

We want to have sufficient PSF variability, so on each visit image we only include one time point for one system. 

And we delete the images after each run to save disk space (maybe keep cutouts or just reconstruction it in the future). 
Just keep the catalog (FITS table).

Note the orientation of the stamp only depends on the boresight angle of the visit. 


---
We first start from a small dataset -- 10x10 stamps on 10'x10'. Each stamp is an independent lens system.
Then each visit corresponds to a time point (epoch) of the lens system. We then extend this to 21x21 stamps on 10'x10', which covers the 1st sample of [simulated lens](https://drive.google.com/drive/u/0/folders/1v4q5IjwjELl8dOaWQQqQ7VGg1YqXC1hd). We currently just test the i-band.
