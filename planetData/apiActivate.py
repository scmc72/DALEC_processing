import os
import requests

just_clip = {
				"name": "just clip",
				"source_type": "scenes",
				"source_type": "scenes",
				"products": [
				{
				  "item_ids": [
				    "20170614_113217_3163208_RapidEye-5"
				  ],
				  "item_type": "REOrthoTile",
				  "product_bundle": "analytic"
				}
				],
				"tools": [
				{
				  "clip": {
				    "aoi": {
				      "type": "Polygon",
				      "coordinates": [
				        [
				          [
				            -163.828125,
				            -44.59046718130883
				          ],
				          [
				            181.7578125,
				            -44.59046718130883
				          ],
				          [
				            181.7578125,
				            78.42019327591201
				          ],
				          [
				            -163.828125,
				            78.42019327591201
				          ],
				          [
				            -163.828125,
				            -44.59046718130883
				          ]
				        ]
				      ]
				    }
				  }
				}
				]
				}


item_id = "20160707_195147_1057916_RapidEye-1"
item_type = "REOrthoTile"
asset_type = "visual"

# setup auth
session = requests.Session()
session.auth = ('PLAK592478481dd34940863ba3b15e2c5a47', '')

# request an item
item = \
  session.get(
    ("https://api.planet.com/data/v1/item-types/" +
    "{}/items/{}/assets/").format(item_type, item_id))

# extract the activation url from the item for the desired asset
item_activation_url = item.json()[asset_type]["_links"]["activate"]

# request activation
response = session.post(item_activation_url)

print(response.status_code)

# send this command in terminal to get status of activation 
# curl -L -H "Authorization: api-key PLAK592478481dd34940863ba3b15e2c5a47" \
#    'https://api.planet.com/data/v1/item-types/REOrthoTile/items/20160707_195147_1057916_RapidEye-1/assets/' \
#    | jq .visual.status

# once active, we can send this to get the download link
# curl -L -H "Authorization: api-key PLAK592478481dd34940863ba3b15e2c5a47" \
#    'https://api.planet.com/data/v1/item-types/REOrthoTile/items/20160707_195147_1057916_RapidEye-1/assets/' \
#    | jq .visual.location