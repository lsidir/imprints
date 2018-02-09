library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity binary_search is
generic(
	COMP_WIDTH : integer := 8;
	POS_WIDTH : integer := 6;     -- 6 bits for position 2^6 = 64 number of buckets
	DATA_WIDTH : integer := 64;   -- 64 bits data width
	LOG_NUM_COMP : integer := 3); -- log(number) of comparators
port(
	clk : in std_logic;
	resetn : in std_logic;

	in_data : in std_logic_vector(DATA_WIDTH-1 downto 0);                      -- input data of DATA_WIDTH bits
	in_data_valid : in std_logic;                                              -- input data is valid
	in_pos : in std_logic_vector(POS_WIDTH-1 downto 0);                        -- input position to add

	in_config : in std_logic_vector(2**LOG_NUM_COMP*DATA_WIDTH-1 downto 0);    -- input configuration (bounds of bins in imprints)
	in_config_valid : in std_logic_vector(COMP_WIDTH-1 downto 0);              -- which input configuration is valid

	out_data : out std_logic_vector(DATA_WIDTH-1 downto 0);                    -- output data is same as input
	out_data_valid : out std_logic;                                            -- output data valid same as input
	out_pos : out std_logic_vector(POS_WIDTH-1 downto 0);                      -- output new position (input position + new result)
	out_config : out std_logic_vector(2**LOG_NUM_COMP*DATA_WIDTH-1 downto 0);   -- input configuration (bounds of bins in imprints)
	out_config_valid : out std_logic_vector(COMP_WIDTH-1 downto 0));             -- which input configuration is valid
end binary_search;

architecture structure of binary_search is

type inter_data_type is array (COMP_WIDTH downto 0) of std_logic_vector(DATA_WIDTH-1 downto 0);
signal inter_data : inter_data_type;
signal inter_data_valid : std_logic_vector(COMP_WIDTH downto 0);
type inter_pos_type is array (COMP_WIDTH downto 0) of std_logic_vector(POS_WIDTH-1 downto 0);
signal inter_pos : inter_pos_type;

type inter_config_type is array (COMP_WIDTH-1 downto 0) of std_logic_vector(2**LOG_NUM_COMP*DATA_WIDTH-1 downto 0);
signal inter_config : inter_config_type;
signal inter_config_valid : std_logic_vector(COMP_WIDTH-1 downto 0);

component comparator
generic (
	POS_WIDTH : integer := 6;     -- 6 bits for position 2^6 = 64 number of buckets
	DATA_WIDTH : integer := 64;   -- 64 bits data width
	LOG_NUM_COMP : integer := 3); -- log(number) of comparators
port (
	clk : in std_logic;           -- clokc
	resetn : in std_logic;        -- reset
	in_data : in std_logic_vector(DATA_WIDTH-1 downto 0);                      -- input data of DATA_WIDTH bits
	in_data_valid : in std_logic;                                              -- input data is valid
	in_pos : in std_logic_vector(POS_WIDTH-1 downto 0);                        -- input position to add
	in_config : in std_logic_vector(2**LOG_NUM_COMP*DATA_WIDTH-1 downto 0);    -- input configuration (bounds of bins in imprints)
	in_config_valid : in std_logic;                                            -- input configuration is valid

	out_data : out std_logic_vector(DATA_WIDTH-1 downto 0);                    -- output data is same as input
	out_data_valid : out std_logic;                                            -- output data valid same as input
	out_pos : out std_logic_vector(POS_WIDTH-1 downto 0));                     -- output new position (input position + new result)
end component;

begin

inter_data(0) <= in_data;
out_data <= inter_data(COMP_WIDTH);
inter_data_valid(0) <= in_data_valid;
out_data_valid <= inter_data_valid(COMP_WIDTH);
inter_pos(0) <= in_pos;
out_pos <= inter_pos(COMP_WIDTH);

GenCOMP: for k in 0 to COMP_WIDTH-1 generate
		comparatorX : comparator generic map (
									POS_WIDTH => 6,
									DATA_WIDTH => 64,
									LOG_NUM_COMP =>3)
								 port map (
									clk => clk,
									resetn => resetn,
									in_data => inter_data(k),
									in_data_valid => inter_data_valid(k),
									in_pos => inter_pos(k),
									in_config => inter_config(k),
									in_config_valid => inter_config_valid(k),
									out_data => inter_data(k+1),
									out_data_valid => inter_data_valid(k+1),
									out_pos => inter_pos(k+1));
end generate GenCOMP;

process(clk)
begin
	if clk'event and clk = '1' then
		out_config_valid <= in_config_valid;
		out_config <= in_config;
		inter_config_valid <= in_config_valid;
		inter_config <= (others => in_config);
	end if; -- clk'event and clk = '1'
end process;

end structure;
