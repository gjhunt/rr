run:
	echo "Making data matrices..."
	Rscript make_data_matrices.R
	echo "Making factors..."
	Rscript make_factors.R
	echo "Estimating..."
	Rscript estimate.R -tag=_no_nid_eln
	Rscript estimate.R
