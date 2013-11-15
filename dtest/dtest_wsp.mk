.PHONY: clean All

All:
	@echo "----------Building project:[ itertest - Debug ]----------"
	@cd "itertest" && $(MAKE) -f  "itertest.mk"
clean:
	@echo "----------Cleaning project:[ itertest - Debug ]----------"
	@cd "itertest" && $(MAKE) -f  "itertest.mk" clean
