library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity imprint_creator is
generic(ADDR_WIDTH : integer := 32;
		DATA_WIDTH : integer := 64);       -- 64 bits data width
port(
	clk : in std_logic;
	resetn : in std_logic;

	read_request : out std_logic;
	read_request_address : out std_logic_vector(ADDR_WIDTH-1 downto 0);
	read_request_almostfull : in std_logic;

	read_response : in std_logic;
	read_response_data : in std_logic_vector(511 downto 0);
	read_response_address : in std_logic_vector(ADDR_WIDTH-1 downto 0);

	write_request : out std_logic;
	write_request_address : out std_logic_vector(ADDR_WIDTH-1 downto 0);
	write_request_data : out std_logic_vector(511 downto 0);
	write_request_almostfull : in std_logic;

	write_response : in std_logic;

	start : in std_logic;
	done : out std_logic;

	-- Parameters
	pos_width : in std_logic(5 downto 0);
end imprint_creator;

architecture structure of imprint_creator is

component binary_search
generic (
port (
end component;

begin


GenCOMP:
	for k in 0 to COMP_HEIGHT-1 generate
end generate GenCOMP;

process(clk)
begin

if clk'event and clk = '1' then
	if resetn = '0' then
		read_request <= '0';
		read_request_address <= (others => '0');
	else -- resetn != '0'

	end if; -- resetn = '0'
end if; -- clk'event and clk = '1'
end process;
end structure;
