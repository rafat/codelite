.PHONY: clean All

All:
	@echo "----------Building project:[ newtonmin - Debug ]----------"
	@cd "newtonmin" && $(MAKE) -f  "newtonmin.mk"
clean:
	@echo "----------Cleaning project:[ newtonmin - Debug ]----------"
	@cd "newtonmin" && $(MAKE) -f  "newtonmin.mk" clean
