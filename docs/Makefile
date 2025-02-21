THEME_VERSION := v0.13.2
THEME := hugo-geekdoc
DESTINATION := target-doc
THEMEDIR := $(DESTINATION)/themes
CURRENTDIR := $(shell pwd)

.PHONY: serve
serve:
	cd $(DESTINATION); hugo serve  --config=$(CURRENTDIR)/config.toml  --contentDir=$(CURRENTDIR)/content/

.PHONY: doc
doc: clean doc-assets doc-build

.PHONY: doc-assets
doc-assets:
	mkdir -p $(THEMEDIR)/$(THEME)/ ; \
	curl -sSL "https://github.com/thegeeklab/$(THEME)/releases/download/${THEME_VERSION}/$(THEME).tar.gz" | tar -xz -C $(THEMEDIR)/$(THEME)/ --strip-components=1

.PHONY: doc-build
doc-build: doc-assets
	cd $(DESTINATION); hugo --minify --config=$(CURRENTDIR)/config.toml  --contentDir=$(CURRENTDIR)/content/

.PHONY: doc-deploy
doc-deploy: clean doc-assets doc-deploy
	cd $(DESTINATION); hugo --minify --config=$(CURRENTDIR)/config.toml  --contentDir=$(CURRENTDIR)/content/ --baseURL https://nokx5.github.io/golden_cpp/

.PHONY: clean
clean:
	rm -rf $(THEMEDIR) && \
	rm -rf $(DESTINATION)/public
