.PHONY: plots analysis

plots:
	echo "Cleaning up..."
	find ./MCF10A/ -name "*.pdf" -exec rm {} \;
	find ./MCF10A/ -name "*.eps" -exec rm {} \;
	echo "Plotting..."
	cd ./MCF10A/plot_scripts && Rscript plot.R
	cd ./MCF10A/plot_scripts && Rscript plot.R -sffx=_no_nid_eln
	rsync -av --exclude='*.R*' /MCF10A/plots /output
	rsync -av --exclude='*.R*' /MCF10A/plots_no_nid_eln /output

analysis:
	echo "Cleaning up..."
	rm -r ./MCF10A/data/ || true
	rm -r ./MCF10A/scaledY/ || true
	rm -r ./MCF10A/plots* || true
	echo "Running processing scripts..."
	cd ./MCF10A/processing_scripts && make
	echo "Cleaning up..."
	find ./MCF10A/ -name "*.pdf" -exec rm {} \;
	find ./MCF10A/ -name "*.eps" -exec rm {} \;
	echo "Plotting..."
	cd ./MCF10A/plot_scripts && Rscript plot.R
	cd ./MCF10A/plot_scripts && Rscript plot.R -sffx=_no_nid_eln
	rsync -av --exclude='*.R*' /MCF10A/plots /output
	rsync -av --exclude='*.R*' /MCF10A/plots_no_nid_eln /output


