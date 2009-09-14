  INCLUDE(BundleUtilities)
  clear_bundle_keys(keys)
# Apparently no plugins are required
#  file(GLOB pluginlist "${QT_PLUGINS_DIR}/imageformats/*.dylib")
  get_bundle_keys(${KNUT_BINARY_DIR}/gui/Knut.app "${pluginlist}" "${KNUT_BINARY_DIR}/Knut.app/Contents/Plugins/" keys)
  fixup_bundle(${KNUT_BINARY_DIR}/Knut.app "" "${KNUT_BINARY_DIR}/Knut.app/Contents/Plugins/")
  foreach(key ${keys})
        execute_process(COMMAND lipo ${${key}_RESOLVED_EMBEDDED_ITEM} -thin ${OSX_ARCH} -output ${${key}_RESOLVED_EMBEDDED_ITEM})
  endforeach(key)
