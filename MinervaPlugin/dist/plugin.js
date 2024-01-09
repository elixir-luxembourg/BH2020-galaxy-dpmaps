(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
var css = ".starter-kit-container {\n  padding: 10px;\n  background-color: #ffffff;\n}\n.disabledbutton {\n  pointer-events: none;\n  opacity: 0.4;\n}\ninput[type=checkbox] .cbcontainer {\n  margin-left: 0;\n}\nhr {\n  border-radius: 5px;\n  background: linear-gradient(to right, #f28988, #6c6a48);\n}\n.custom_checkbox {\n  margin: 2px 0px 0px 0px;\n  text-align: center !important;\n  color: #6d6d6d;\n  cursor: pointer;\n}\n.custom_checkbox,\ninput[type=checkbox] {\n  vertical-align: middle;\n}\n.custom_checkbox_label {\n  font-family: lato;\n  color: #333;\n}\n.custom_select_label {\n  padding-top: 4px !important;\n  padding: 0;\n  text-align: right;\n  margin-left: 5px;\n}\n.inputfile {\n  font-size: 1.0em;\n  font-weight: 500;\n  color: black;\n  display: inline-block;\n}\n.btn-download {\n  background-color: DodgerBlue;\n  border: none;\n  color: white;\n  padding: 12px 30px;\n  cursor: pointer;\n  font-size: 20px;\n}\n.btn-download:hover {\n  background-color: RoyalBlue;\n}\n.btn-group {\n  width: 100%;\n}\n.clickElementinTable {\n  background-color: transparent;\n  display: inline-block;\n  white-space: nowrap;\n}\n.clickElementinTable body {\n  text-align: center;\n}\n"; (require("browserify-css").createStyle(css, { "href": "BH2020_Plugin\\__tmp_minerva_plugin\\css\\styles.css" }, { "insertAt": "bottom" })); module.exports = css;
},{"browserify-css":3}],2:[function(require,module,exports){
require('../css/styles.css');

const pluginName = 'Galaxy';
const pluginVersion = '0.9.0';
const minervaProxyServer = 'https://minerva-dev.lcsb.uni.lu/minerva-proxy/';
let GalaxyData = [];
let GalaxyTable = null;
let maxFC = 0;
let Header = true;
let CSV = false;
const globals = {
  selected: [],
  allBioEntities: [],
  pickedRandomly: undefined,
  container: undefined,
  Overlays: ['InDegree', 'OutDegree', 'Betweenness', 'Closeness'],
  allSpeciesNames: {},
  allSpecies: {},
  allReactions: {},
  allHGNC: {},
  allEntrez: {},
  allUniprot: {},
  mapping: "",
  count: 0,
  downloadtext: '',
  checkbox_checked: undefined,
  modelid: -1,
  progress: 0,
  firstelement: undefined,
  secondelement: undefined,
  phenotype: [],
  element: undefined
};
let $ = window.$;

if ($ === undefined && minerva.$ !== undefined) {
  $ = minerva.$;
} // ******************************************************************************
// ********************* PLUGIN REGISTRATION WITH MINERVA *********************
// ******************************************************************************


let minervaProxy;
let pluginContainer;
let pluginContainerId;
let minervaVersion;

const register = function (_minerva) {
  console.log('registering ' + pluginName + ' plugin');
  $(".tab-content").css('position', 'relative');
  minervaProxy = _minerva;
  pluginContainer = $(minervaProxy.element);
  pluginContainerId = pluginContainer.attr('id');

  if (!pluginContainerId) {
    //the structure of plugin was changed at some point and additional div was added which is the container but does not have any properties (id or height)
    pluginContainerId = pluginContainer.parent().attr('id');
  }

  console.log('minerva object ', minervaProxy);
  console.log('project id: ', minervaProxy.project.data.getProjectId());
  console.log('model id: ', minervaProxy.project.data.getModels()[0].modelId);
  return minerva.ServerConnector.getConfiguration().then(function (conf) {
    minervaVersion = parseFloat(conf.getVersion().split('.').slice(0, 2).join('.'));
    console.log('minerva version: ', minervaVersion);
    initPlugin();
  });
};

const unregister = function () {
  console.log('unregistering ' + pluginName + ' plugin');
  unregisterListeners();
  return deHighlightAll();
};

const getName = function () {
  return pluginName;
};

const getVersion = function () {
  return pluginVersion;
};
/**
 * Function provided by Minerva to register the plugin
 */


minervaDefine(function () {
  return {
    register: register,
    unregister: unregister,
    getName: getName,
    getVersion: getVersion,
    minWidth: 400,
    defaultWidth: 500
  };
});

function initPlugin() {
  registerListeners();
  initMainPageStructure();
}

function registerListeners() {
  minervaProxy.project.map.addListener({
    dbOverlayName: "search",
    type: "onSearch",
    callback: searchListener
  });
}

function unregisterListeners() {
  minervaProxy.project.map.removeAllListeners();
} // ****************************************************************************
// ********************* MINERVA INTERACTION*********************
// ****************************************************************************


function deHighlightAll() {
  return minervaProxy.project.map.getHighlightedBioEntities().then(highlighted => minervaProxy.project.map.hideBioEntity(highlighted));
} // ****************************************************************************
// ********************* PLUGIN STRUCTURE AND INTERACTION*********************
// ****************************************************************************


function initMainPageStructure() {
  globals.container = $('<div class="' + pluginName + '-container"></div>').appendTo(pluginContainer);
  $(`<div id="gal_stat_spinner" class="mt-5">
        
    </div>`).appendTo(globals.container);
  var url = new URL(window.location.href);
  var query = url.searchParams.get("datasource");

  if (query == null) {
    $("#gal_stat_spinner").html(`
            <div class="alert alert-danger ml-2 mr-2" role="alert">
                <span><i class="fas fa-exclamation-triangle"></i></span>
                <span class="sr-only">Error:</span>
                Please supply a 'datasource' parameter within the URL.
            </div>    
        `);
    return;
  }

  try {
    new URL(query);
  } catch (_) {
    $("#gal_stat_spinner").html(`
            <div class="alert alert-danger ml-2 mr-2" role="alert">
                <span><i class="fas fa-exclamation-triangle"></i></span>
                <span class="sr-only">Error:</span>
                The data source is not a valid URL.
            </div>    
        `);
    return;
  }

  $("#gal_stat_spinner").html(`
        <div class="d-flex justify-content-center">
            <div class="spinner-border" role="status">
                <span class="sr-only"></span>
            </div>
        </div>
        <div class="d-flex justify-content-center mt-2">
            <span id="gal_loading_text">Fetching Data ...</span>
        </div>
    `);

  try {
    Header = url.searchParams.get("header") != "false";
  } catch (_) {}

  try {
    CSV = url.searchParams.get("datatype") == "csv";
  } catch (_) {}

  fetchGalaxyQuery(query).then(raw_galaxy_data => {
    let mapping_dict = {};

    switch (globals.mapping) {
      case "identifier_hgnc_symbol":
        mapping_dict = globals.allHGNC;
        break;

      case "identifier_entrez":
        mapping_dict = globals.allEntrez;
        break;

      case "identifier_uniprot":
        mapping_dict = globals.allUniprot;
        break;

      default:
        mapping_dict = globals.allHGNC;
        break;
    }

    $("#gal_loading_text").html("Reading Map Elements ...");
    minervaProxy.project.data.getAllBioEntities().then(function (bioEntities) {
      $("#gal_loading_text").html("Generating Table ...");
      globals.allBioEntities = bioEntities;

      for (var e of globals.allBioEntities) {
        if (e.constructor.name === 'Alias') {
          var ename = e.getName().toLowerCase();

          if (!globals.allSpeciesNames.hasOwnProperty(ename)) {
            globals.allSpeciesNames[ename] = [];
          }

          globals.allSpeciesNames[ename].push(e);
          globals.allSpecies[e.id] = e;

          for (let reference of e.references) {
            var reference_id = reference._resource.toLowerCase();

            switch (reference._type) {
              case "ENTREZ":
                if (!globals.allEntrez.hasOwnProperty(reference_id)) {
                  globals.allEntrez[reference_id] = [];
                }

                globals.allEntrez[reference_id].push(e);
                break;

              case "HGNC_SYMBOL":
                if (!globals.allHGNC.hasOwnProperty(reference_id)) {
                  globals.allHGNC[reference_id] = [];
                }

                globals.allHGNC[reference_id].push(e);
                break;

              case "UNIPROT":
                if (!globals.allUniprot.hasOwnProperty(reference_id)) {
                  globals.allUniprot[reference_id] = [];
                }

                globals.allUniprot[reference_id].push(e);
                break;

              default:
                break;
            }
          }
        } else {
          globals.allReactions[e.id] = e;
        }
      }

      ;

      for (var [id, data] of Object.entries(raw_galaxy_data)) {
        if (mapping_dict.hasOwnProperty(id)) {
          var mapped_minerva_elements = mapping_dict[id];

          if (mapped_minerva_elements.length > 0) {
            data["name"] = mapped_minerva_elements[0].getName();
            data["minerva_elements"] = mapped_minerva_elements;
            GalaxyData.push(data);
          }

          if (Math.abs(data.fc) > maxFC) {
            maxFC = Math.abs(data.fc);
          }
        }
      }

      document.getElementById("gal_stat_spinner").remove();
      globals.container.append(`

                <div class="cbcontainer mt-4 mb-2">
                    <input type="checkbox" class="custom_checkbox" id="gal_checkbox_fc_enable">
                    <label class="custom_checkbox" for="gal_checkbox_fc_enable">Use FC Threshold:</label>
                </div>
                <div class="row mb-4 mt-2 ml-2 disabledbutton" id="gal_fc-container">
                    <div class="col-auto custom_select_label" style="padding:0; width: 30%; text-align: right; ">
                        <span style="margin: 0; display: inline-block; vertical-align: middle; line-height: normal;">Absolute FC Threshold:</span>
                    </div>
                    <div class="col">
                        <input type="text" class="textfield" value="1.00" id="gal_fcthreshold"/>
                    </div>
                </div>

                
                <div class="cbcontainer mt-4 mb-2">
                    <input type="checkbox" class="custom_checkbox" id="gal_checkbox_pvalue_enable" checked>
                    <label class="custom_checkbox" for="gal_checkbox_pvalue_enable">Use p-value Threshold:</label>
                </div>

                <div class="row mt-2 mb-4 ml-2" id="gal_pvalue-container">
                    <div class="col-auto custom_select_label" style="padding:0; width: 30%; text-align: right; ">
                        <span style="margin: 0; display: inline-block; vertical-align: middle; line-height: normal;">p-value Threshold:</span>
                    </div>
                    <div class="col">
                        <input type="text" class="textfield" value="0.05" id="gal_pvaluethreshold"/>
                    </div>
                    <div class="col">
                        <div class="cbcontainer">
                            <input type="checkbox" class="custom_checkbox" id="gal_checkbox_adjusted" checked>
                            <label class="custom_checkbox" for="checkbox_adjusted">Adjusted p-value?</label>
                        </div>
                    </div>
                </div>


                <div class="mt-4 mb-4">
                    <table style="width:100%" class="table nowrap table-sm" id="gal_galaxy_table" cellspacing="0">
                        <thead>
                            <tr>
                                <th style="vertical-align: middle;"></th>
                                <th style="vertical-align: middle;">Element</th>
                                <th style="vertical-align: middle;">FC-Value</th>
                                <th style="vertical-align: middle;">p-value</th>
                                <th style="vertical-align: middle;">adj. p-value</th>
                            </tr>
                        </thead>
                    </table>
                </div>
                <button type="button" id="gal_btn_selectall" class="btn btn-primary btn-block mt-4 mb-2">Select All</button>
                <button type="button" id="gal_btn_reset" class="btn btn-primary btn-block mt-2">Reset</button>
            `);
      $("#gal_pvaluethreshold").on("input", fillTable);
      $("#gal_fcthreshold").on("input", fillTable);
      $("#gal_checkbox_adjusted").on("change", fillTable);
      $("#gal_checkbox_pvalue_enable").on("change", function () {
        if ($(this).prop('checked') === true) {
          $("#gal_pvalue-container").removeClass("disabledbutton");
        } else {
          $("#gal_pvalue-container").addClass("disabledbutton");
        }

        fillTable();
      });
      $("#gal_checkbox_fc_enable").on("change", function () {
        if ($(this).prop('checked') === true) {
          $("#gal_fc-container").removeClass("disabledbutton");
        } else {
          $("#gal_fc-container").addClass("disabledbutton");
        }

        fillTable();
      });
      $('#gal_btn_selectall').on('click', async function () {
        var _text = await disablebutton("gal_btn_selectall");

        var elements = [];
        dehighlightall().then(r => {
          GalaxyTable.rows().every(function () {
            var row = this.nodes().to$();
            row.find('.gal_clickCBinTable').prop('checked', true);
            elements.push(row.find('.gal_clickCBinTable').attr("data"));
          });
          highlightElements(elements).finally(r => {
            enablebtn("gal_btn_selectall", _text);
          });
        });
      });
      $("#target").keypress(function () {
        console.log("Handler for .keypress() called.");
      });
      $("#gal_fcthreshold").bind("keypress", function (evt) {
        evt = evt ? evt : window.event;
        var charCode = evt.which ? evt.which : evt.keyCode;

        if (charCode > 31 && (charCode < 48 || charCode > 57) && charCode != 44 && charCode != 46) {
          return false;
        }

        return true;
      });
      $("#gal_pvaluethreshold").bind("keypress", function (evt) {
        evt = evt ? evt : window.event;
        var charCode = evt.which ? evt.which : evt.keyCode;

        if (charCode > 31 && (charCode < 48 || charCode > 57) && charCode != 44 && charCode != 46) {
          return false;
        }

        return true;
      });
      $('#gal_btn_reset').on('click', async function () {
        var _text = await disablebutton("gal_btn_reset");

        dehighlightall().then(r => {
          GalaxyTable.rows().every(function () {
            var row = this.nodes().to$();
            row.find('.gal_clickCBinTable').prop('checked', false);
          });
          enablebtn("gal_btn_reset", _text);
        });
      });
      GalaxyTable = $('#gal_galaxy_table').DataTable({
        "order": [[1, "asc"]],
        scrollX: true,
        autoWidth: true,
        columns: [{
          "width": "10%"
        }, null, {
          "width": "20%"
        }, {
          "width": "20%"
        }, {
          "width": "20%"
        }],
        columnDefs: [{
          targets: 0,
          className: 'dt-center'
        }, {
          targets: 1,
          className: 'dt-center'
        }, {
          targets: 2,
          className: 'dt-center'
        }, {
          targets: 3,
          className: 'dt-center'
        }, {
          targets: 4,
          className: 'dt-center'
        }]
      }).columns.adjust().draw();
      fillTable();
    });
  }).catch(e => {
    $("#gal_stat_spinner").html(`
            <div class="alert alert-danger ml-2 mr-2" role="alert">
                <span><i class="fas fa-exclamation-triangle"></i></span>
                <span class="sr-only">Error:</span>
                Could not fetch data from the supplied URL.
            </div>    
        `);
  });
}

async function disablebutton(id) {
  var promise = new Promise(function (resolve, reject) {
    setTimeout(() => {
      pluginContainer.addClass("disabledbutton");
      var $btn = $('#' + id);
      let text = $btn.html();
      $btn.html('<span class="loadingspinner spinner-border spinner-border-sm"></span>');
      resolve(text);
    }, 0);
  });
  return promise;
}

async function enablebtn(id, text) {
  return new Promise(resolve => {
    setTimeout(() => {
      pluginContainer.removeClass("disabledbutton");
      var $btn = $('#' + id);
      $btn.html(text);
      resolve('');
    }, 0);
  });
}

async function fillTable() {
  let pvalue_threshold = parseFloat($("#gal_pvaluethreshold").val().replace(',', '.'));

  if (isNaN(pvalue_threshold)) {
    pvalue_threshold = 1;
  }

  let fc_threshold = parseFloat($("#gal_fcthreshold").val().replace(',', '.'));

  if (isNaN(pvalue_threshold)) {
    fc_threshold = 0;
  }

  GalaxyTable.clear();

  for (var [id, data] of Object.entries(GalaxyData)) {
    if (isNaN(data.fc) == false && Math.abs(data.fc) < fc_threshold && $("#gal_checkbox_fc_enable")[0].checked) continue;

    if ($("#gal_checkbox_pvalue_enable")[0].checked) {
      if ($('#gal_checkbox_adjusted')[0].checked) {
        if (isNaN(data.pvalue_adj) == false && data.pvalue_adj > pvalue_threshold) continue;
      } else {
        if (isNaN(data.pvalue) == false && data.pvalue > pvalue_threshold) continue;
      }
    }

    var result_row = [`<input type="checkbox" class="gal_clickCBinTable" data="${id}">`, `<a href="#gal_" data="${id}" class="gal_elementlink">${data.name}</a>`, '<span data-toggle="tooltip" title="' + data.fc + '">' + expo(data.fc, 3, 2) + '</span>', '<span data-toggle="tooltip" title="' + data.pvalue + '">' + expo(data.pvalue, 4, 2) + '</span>', '<span data-toggle="tooltip" title="' + data.pvalue_adj + '">' + expo(data.pvalue_adj, 4, 2) + '</span>'];
    GalaxyTable.row.add(result_row);
  }

  GalaxyTable.columns.adjust().draw();
}

$('#gal_galaxy_table').on('draw.dt', function () {
  $('[data-toggle="tooltip"]').tooltip();
});

function fetchGalaxyQuery(query) {
  return new Promise((resolve, reject) => {
    var client = new XMLHttpRequest();
    let mapping = "name";
    client.open('GET', query);

    client.onerror = function () {
      alert("The data does'nt have the right CORS headers, please ask your admin to fix it.");
      reject();
    };

    client.onreadystatechange = function () {
      if (this.readyState == 4) {
        if (this.status == 200) {
          var output = {};
          var raw_response = client.responseText;
          var firstline = true;

          for (var line of raw_response.split("\n")) {
            if (firstline) {
              firstline = false;

              if (Header) {
                globals.mapping = line.split(CSV ? "," : "\t")[0];
                continue;
              }
            }

            var entries = line.split(CSV ? "," : "\t");

            if (entries.length < 2) {
              continue;
            }

            var fc = parseFloat(entries[1]);

            if (isNaN(fc)) {
              fc = "NaN";
            }

            var pvalue = entries.length > 2 ? parseFloat(entries[2]) : "N/A";
            var adj_pvalue = entries.length > 3 ? parseFloat(entries[3]) : "N/A";
            output[entries[0].toLowerCase()] = {
              "id": entries[0],
              "fc": fc,
              "pvalue": pvalue,
              "pvalue_adj": adj_pvalue
            };
          }

          resolve(output);
        } else {
          console.log(client.responseText);
          reject();
        }
      }
    };

    client.send();
  });
}

function dehighlightElement(element) {
  var minerva_elements = GalaxyData[element].minerva_elements;
  minervaProxy.project.map.getHighlightedBioEntities().then(highlighted => {
    minervaProxy.project.map.hideBioEntity(minerva_elements.filter(e => highlighted.includes(e)));
  });
}

function dehighlightall() {
  return new Promise((resolve, reject) => {
    minervaProxy.project.map.getHighlightedBioEntities().then(highlighted => {
      minervaProxy.project.map.hideBioEntity(highlighted).finally(r => {
        resolve();
      });
    });
  });
}

function highlightElements(node_ids) {
  return new Promise((resolve, reject) => {
    const highlightDefs = [];

    for (var node_id of node_ids) {
      var e = GalaxyData[node_id];
      var _value = 0;
      if (!isNaN(e.fc)) _value = maxFC != 0 ? e.fc / maxFC : 0;
      var hex = rgbToHex((1 - Math.abs(_value)) * 255);
      if (_value > 0) hex = '#ff' + hex + hex;else if (_value < 0) hex = '#' + hex + hex + 'ff';else hex = '#ffffff';

      for (var minerva_element of e.minerva_elements) {
        highlightDefs.push({
          element: {
            id: minerva_element.id,
            modelId: minerva_element.getModelId(),
            type: minerva_element.constructor.name.toUpperCase()
          },
          type: "SURFACE",
          options: {
            color: hex,
            opacity: 0.6
          }
        });
      }
    }

    minervaProxy.project.map.showBioEntity(highlightDefs).then(r => {
      resolve();
    });
  });
}

function rgbToHex(rgb) {
  var hex = Number(Math.round(rgb)).toString(16);

  if (hex.length < 2) {
    hex = "0" + hex;
  }

  return hex;
}

;

function highlightSelected(reset = true, _name) {
  minervaProxy.project.map.getHighlightedBioEntities().then(highlighted => {
    minervaProxy.project.map.hideBioEntity(reset ? highlighted : []).then(r => {
      const highlightDefs = [];

      for (var _name in globals.selected_elements) {
        for (let _id of globals.allSpeciesNames[_name]) {
          var e = globals.allSpecies[_id];
          highlightDefs.push({
            element: {
              id: e.id,
              modelId: e.getModelId(),
              type: e.constructor.name.toUpperCase()
            },
            type: "SURFACE",
            options: {
              color: '#gal_FF1593',
              opacity: 0.5
            }
          });
        }
      }
      /*
      if (pickedRandomly) {
          if (globals.pickedRandomly) {
              highlightDefs.push({
                  element: {
                      id: globals.pickedRandomly.id,
                      modelId: globals.pickedRandomly.getModelId(),
                      type: globals.pickedRandomly.constructor.name.toUpperCase()
                  },
                  type: "SURFACE",
                  options: {
                      color: '#gal_00FF00',
                      opacity: 0.2
                  }
              });
          }
      } else {
          globals.selected.forEach(e => {
              if (e.constructor.name === 'Alias') {
                  highlightDefs.push({
                      element: {
                          id: e.id,
                          modelId: e.getModelId(),
                          type: "ALIAS"
                      },
                      type: "ICON"
                  });
              }
          });
      }
      */


      minervaProxy.project.map.showBioEntity(highlightDefs);
    });
  });
}

function focusOnSelected() {
  function focus(entity) {
    if (entity.constructor.name === 'Alias') {
      minervaProxy.project.map.fitBounds({
        modelId: entity.getModelId(),
        x1: entity.getX(),
        y1: entity.getY(),
        x2: entity.getX() + entity.getWidth(),
        y2: entity.getY() + entity.getHeight()
      });
    } else {
      minervaProxy.project.map.fitBounds({
        modelId: entity.getModelId(),
        x1: entity.getCenter().x,
        y1: entity.getCenter().y,
        x2: entity.getCenter().x,
        y2: entity.getCenter().y
      });
    }
  }

  if (globals.selected.length > 0) {
    minervaProxy.project.map.openMap({
      id: globals.selected[0].getModelId()
    });
    focus(globals.selected[0]);
  }
}

$(document).on('click', '.gal_elementlink', function () {
  selectElementonMap($(this).attr('data'), false);
});
$(document).on('change', '.gal_clickCBinTable', function () {
  if ($(this).prop('checked') === true) {
    highlightElements([$(this).attr('data')]);
  } else {
    dehighlightElement($(this).attr('data'));
  }
});

function selectElementonMap(element) {
  var minerva_elements = GalaxyData[element].minerva_elements;
  globals.selected = [];

  if (minerva_elements.length > 0) {
    globals.selected.push(minerva_elements[0]);
    focusOnSelected();
  }
}

function searchListener(entites) {
  globals.selected = entites[0];
  let str = '';

  if (globals.selected.length > 0) {
    globals.selected.forEach(e => {
      if (e.constructor.name === 'Alias') str += `<div>${e.getName()} - ${e.getElementId()} - ${e._type}</div>`;
    });
  }

  pluginContainer.find('.panel-test .panel-body').html(str);
}

function expo(x, f = 3, e = 3) {
  let _round = Math.floor(x * Math.pow(10, f - 1)) / Math.pow(10, f - 1);

  if (_round == 0) return x.toExponential(e);else return Math.round(x * Math.pow(10, f)) / Math.pow(10, f);
}
},{"../css/styles.css":1}],3:[function(require,module,exports){
'use strict';
// For more information about browser field, check out the browser field at https://github.com/substack/browserify-handbook#browser-field.

var styleElementsInsertedAtTop = [];

var insertStyleElement = function(styleElement, options) {
    var head = document.head || document.getElementsByTagName('head')[0];
    var lastStyleElementInsertedAtTop = styleElementsInsertedAtTop[styleElementsInsertedAtTop.length - 1];

    options = options || {};
    options.insertAt = options.insertAt || 'bottom';

    if (options.insertAt === 'top') {
        if (!lastStyleElementInsertedAtTop) {
            head.insertBefore(styleElement, head.firstChild);
        } else if (lastStyleElementInsertedAtTop.nextSibling) {
            head.insertBefore(styleElement, lastStyleElementInsertedAtTop.nextSibling);
        } else {
            head.appendChild(styleElement);
        }
        styleElementsInsertedAtTop.push(styleElement);
    } else if (options.insertAt === 'bottom') {
        head.appendChild(styleElement);
    } else {
        throw new Error('Invalid value for parameter \'insertAt\'. Must be \'top\' or \'bottom\'.');
    }
};

module.exports = {
    // Create a <link> tag with optional data attributes
    createLink: function(href, attributes) {
        var head = document.head || document.getElementsByTagName('head')[0];
        var link = document.createElement('link');

        link.href = href;
        link.rel = 'stylesheet';

        for (var key in attributes) {
            if ( ! attributes.hasOwnProperty(key)) {
                continue;
            }
            var value = attributes[key];
            link.setAttribute('data-' + key, value);
        }

        head.appendChild(link);
    },
    // Create a <style> tag with optional data attributes
    createStyle: function(cssText, attributes, extraOptions) {
        extraOptions = extraOptions || {};

        var style = document.createElement('style');
        style.type = 'text/css';

        for (var key in attributes) {
            if ( ! attributes.hasOwnProperty(key)) {
                continue;
            }
            var value = attributes[key];
            style.setAttribute('data-' + key, value);
        }

        if (style.sheet) { // for jsdom and IE9+
            style.innerHTML = cssText;
            style.sheet.cssText = cssText;
            insertStyleElement(style, { insertAt: extraOptions.insertAt });
        } else if (style.styleSheet) { // for IE8 and below
            insertStyleElement(style, { insertAt: extraOptions.insertAt });
            style.styleSheet.cssText = cssText;
        } else { // for Chrome, Firefox, and Safari
            style.appendChild(document.createTextNode(cssText));
            insertStyleElement(style, { insertAt: extraOptions.insertAt });
        }
    }
};

},{}]},{},[2]);
