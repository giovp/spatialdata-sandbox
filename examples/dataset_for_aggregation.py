##
import importlib.util
import sys

def import_from_path(module_name, file_path):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


my_module = import_from_path("my_module", "../spatialdata/tests/conftest.py")

# Now you can use my_module just like any other module
f = getattr(my_module, '_make_sdata_for_testing_querying_and_aggretation')
sdata = f()
print(sdata)

##
# to visualize the cases considered in the test, much more immediate than reading them as text as done above
import matplotlib.pyplot as plt
import spatialdata_plot

# otherwise pre-commits will remove the "unused" import
_ = spatialdata_plot

ax = plt.gca()
sdata.pl.render_shapes(element="values_polygons", color="categorical_in_obs").pl.show(ax=ax)
sdata.pl.render_shapes(element="values_circles", color="categorical_in_obs").pl.show(ax=ax)
sdata.pl.render_shapes(element="by_polygons", na_color=(1.0, 0.7, 0.7, 0.3)).pl.show(ax=ax)
sdata.pl.render_shapes(element="by_circles", na_color=(1.0, 0.7, 0.7, 0.3)).pl.show(ax=ax)
sdata.pl.render_points(color="categorical_in_ddf", size=10.0, palette="tab10").pl.show(ax=ax)
plt.show()
